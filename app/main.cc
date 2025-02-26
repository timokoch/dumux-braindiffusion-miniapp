// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>
#include <type_traits>

#include <chrono>
#include <thread>

#include <dune/alugrid/grid.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>

#include <dumux/io/gridwriter.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/io/container.hh>
#include <dumux/io/chrono.hh>
#include <dumux/io/gnuplotinterface.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/pdesolver.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/nonlinear/leastsquares.hh>

#include "model.hh"

namespace Dumux {

//////////////////////////////////////////////////////////////
//
// We solve the diffusion equation (see model.hh):
//
//     ∂c/∂t - ∇·(D∇c) = q,
//
// which means on the boundary we have to specify the flux:
//
//     -D∇c·n = k (c_b/φ_b - c/φ),
//
// where c_b is the boundary concentration, φ is the porosity,
// and k is the surface conductance. By dividing by φ, we get
// the fluid concentration c_f = c/φ from the
// total concentration c.
//
// To obtain the boundary concentration c_b at time t, two
// strategies are implemented. The parameter
// "Problem.UseCurveFit" in the parameter file determines
// which strategy is used:
//
// 1. Linear interpolation between given data points (false)
// 2. Biexponential curve fit (true)
//
// The source term q is implemented in the `source` method of
// the problem class. Here we have a simple transfer to blood
// compartment:
//
//     q = c kB,
//
// where kB is the transfer coefficient to the blood compartment.
//
// The diffusion tensor D can be isotropic or anisotropic.
// The parameter "Problem.UseIsotropicMeanDiffusion" in the
// parameter file determines which strategy is used. In the
// isotropic case, the diffusion tensor is diagonal with the
// same value in all directions. The mean value from the
// data set is used. For the anisotropic case, the full
// diffusion tensor is used.
//
// The simulation can be run in parallel with MPI.
//
///////////////////////////////////////////////////////////////
template<class TypeTag>
class DiffusionTestProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::LocalView::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Tensor = Dune::FieldMatrix<double, 3, 3>;
    enum class DiffusionType { DataIsotropic, DataTensor, LiteratureIsotropic };
public:
    DiffusionTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // model parametersr read from parameter file
        k_ = getParam<Scalar>("Problem.SurfaceConductance");
        transferCoeffToBlood_ = getParam<Scalar>("Problem.TransferCoeffToBlood");
        porosityBrain_ = getParam<Scalar>("Problem.PorosityBrain", 0.2);
        porositySAS_ = getParam<Scalar>("Problem.PorositySAS", 1.0);

        // the time points at which the boundary data is given
        times_ = {0, 15716, 86001, 172280, 251641};

        // whether to use a biexpontial curve fit or just linear interpolation in between the given data points
        useCurveFit_ = getParam<bool>("Problem.UseCurveFit", false);

        // whether to use isotropic mean diffusion or the full potentially anisotropic diffusion tensor
        const auto diffType = getParam<std::string>("Problem.DiffusionType", "DataTensor");
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
        const auto dofIndex = scv.dofIndex();
        const auto& volVars = elemVolVars[scv];
        const auto insideConcentration = volVars.priVar(0);
        const auto boundaryConcentration = concentrationBoundaryData(dofIndex);

        // to convert total concentration to fluid concentration, divide by porosity
        return {
            -k_*(boundaryConcentration/porositySAS_ - insideConcentration/porosityBrain_)
        };
    }

    // transfer to blood
    template<class ElementVolumeVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        return {
            -(elemVolVars[scv].priVar(0)/porosityBrain_)*transferCoeffToBlood_
        };
    }

    Tensor diffusionCoefficient(const Element& element, const SubControlVolume& scv) const
    { return diffusionCoefficient_[scv.elementIndex()]; }

    void setTime(const Scalar t) { time_ = t; }

    template<class GridData>
    void setGridData(const GridData& gridData) { setGridData_(gridData); }

    void updateBoundaryDataForOutput(std::vector<Scalar>& cBoundary, const Scalar t)
    {
        const auto& gridView = this->gridGeometry().gridView();
        for (const auto& vertex : vertices(gridView))
        {
            const auto vIdx = this->gridGeometry().dofMapper().index(vertex);
            if (this->gridGeometry().dofOnBoundary(vIdx))
                cBoundary[vIdx] = concentrationBoundaryData(vIdx);
        }
    }

    const Dune::BlockVector<double>& checkPoints() const
    { return times_; }

private:
    Scalar concentrationBoundaryData(int dofIndex) const
    {
        const auto& p = concentrationBoundaryDataParams_[dofIndex];
        return concentrationBoundaryDataFunction_(time_, p);
    }

    Scalar concentrationBoundaryDataFunction_(const Scalar time, const Dune::BlockVector<Scalar>& p) const
    {
        // check if any parameters are not finite
        for (int i = 0; i < p.size(); ++i)
            if (!std::isfinite(p[i]))
                DUNE_THROW(Dune::Exception, "Invalid parameter values: " << p);

        if (useCurveFit_)
        {
            const auto days = time/86400.0;
            return p[0]*days*std::exp(-p[1]*days);
        }
        else
            return Dumux::interpolate<InterpolationPolicy::LinearTable>(time, times_, p);
    }

private:
    Scalar k_;
    Scalar transferCoeffToBlood_;
    Scalar porosityBrain_, porositySAS_;
    std::vector<Tensor> diffusionCoefficient_;
    std::vector<Dune::BlockVector<Scalar>> concentrationBoundaryDataParams_;
    std::vector<Dune::BlockVector<Scalar>> concentrationElementDataParams_;
    Scalar time_;
    bool useCurveFit_;
    DiffusionType diffusionType_;

    Dune::BlockVector<double> times_;

    template<class GridData>
    void setGridData_(const GridData& gridData)
    {
        const auto& gridView = this->gridGeometry().gridView();
        diffusionCoefficient_.resize(gridView.size(0), Tensor(0.0));

        for (const auto& element : elements(gridView))
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            if (diffusionType_ == DiffusionType::DataTensor)
            {
                static const auto waterToGadobutrolRatio = getParam<Scalar>("Problem.DiffusionCoefficientWaterToGadobutrolRatio", 3.5e-10/3e-9); // mm^2/s
                const auto data = gridData.template getArrayParameter<double, 9>(element, "dt");
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        diffusionCoefficient_[eIdx][i][j] = data[i*3 + j] * waterToGadobutrolRatio;
            }
            else if (diffusionType_ == DiffusionType::DataIsotropic)
            {
                static const auto waterToGadobutrolRatio = getParam<Scalar>("Problem.DiffusionCoefficientWaterToGadobutrolRatio", 3.5e-10/3e-9); // mm^2/s
                const auto data = gridData.getParameter(element, "md");
                for (int i = 0; i < 3; ++i)
                    diffusionCoefficient_[eIdx][i][i] = data * waterToGadobutrolRatio;
            }
            else
            {
                const auto marker = gridData.getParameter(element, "subdomains");
                static const auto grayMatterEffDiffusionCoefficient = getParam<Scalar>("Problem.GrayMatterEffectiveDiffusionCoefficient", 1.4e-4); // mm^2/s
                static const auto whiteMatterEffDiffusionCoefficient = getParam<Scalar>("Problem.WhiteMatterEffectiveDiffusionCoefficient", 1.1e-4); // mm^2/s
                for (int i = 0; i < 3; ++i)
                    diffusionCoefficient_[eIdx][i][i] = marker == 1 ? grayMatterEffDiffusionCoefficient : whiteMatterEffDiffusionCoefficient;
            }
        }

        Dune::BlockVector<Scalar> boundaryData(times_.size());
        concentrationBoundaryDataParams_.resize(this->gridGeometry().numDofs());
        std::size_t converged = 0;
        std::size_t numBoundaryVertices = 0;
        for (const auto& vertex : vertices(gridView))
        {
            const auto vIdx = this->gridGeometry().dofMapper().index(vertex);
            if (!this->gridGeometry().dofOnBoundary(vIdx))
                continue;

            ++numBoundaryVertices;

            for (int i = 0; i < times_.size(); ++i)
            {
                boundaryData[i]
                    = gridData.getParameter(vertex, Fmt::format("boundary-concentration-{:d}", int(times_[i])));

                if (!std::isfinite(boundaryData[i]))
                    DUNE_THROW(Dune::Exception, "Invalid boundary data at vertex " << vIdx << " at time " << times_[i]);
            }

            if (useCurveFit_)
            {
                Dune::BlockVector<double> p({0.2, 1.0});
                Dune::BlockVector<double> res(times_.size());
                const auto residual = [&](const auto& params)
                {
                    for (size_t i = 0; i < times_.size(); ++i)
                    {
                        res[i] = this->concentrationBoundaryDataFunction_(times_[i], params) - boundaryData[i];

                        // enforce some bounds on the parameters
                        if (params[0] < 0.0 || params[0] > 10.0)
                            res[i] = params[0]*1e4;

                        if (params[1] < 0.0 || params[1] > 5.0)
                            res[i] = params[1]*1e4;

                        if (!std::isfinite(res[i]))
                            DUNE_THROW(Dune::Exception, "Invalid residual value at vertex " << vIdx << " at time "
                                        << times_[i] << ": " << res[i] << " for parameters " << params);
                    }
                    return res;
                };

                {

                    std::ostringstream local;
                    auto cout_buff = std::cout.rdbuf();
                    std::cout.rdbuf(local.rdbuf());

                    const auto c = Optimization::makeNonlinearLeastSquaresSolver(residual, p, times_.size())->apply(p);

                    std::cout.rdbuf(cout_buff);

                    if (c) ++converged;
                }

                concentrationBoundaryDataParams_[vIdx] = p;
            }
            else
                concentrationBoundaryDataParams_[vIdx] = boundaryData;
        }

        if (useCurveFit_)
            std::cout << "Curve fit converged for " << converged
                      << " / " << numBoundaryVertices << std::endl;
    }
};

} // end namespace Dumux

namespace Dumux::Properties::TTag {

// compile-time configuration of simulation
struct DiffusionTest
{
    using InheritsFrom = std::tuple<DiffusionModel, BoxModel>;

    using Scalar = double;
    using Grid = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;

    template<class TypeTag>
    using Problem = DiffusionTestProblem<TypeTag>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

} // end namespace Dumux::Properties::TTag


//////////////////////////////////////////////////////////////
//
// Main function
//
// The main function initializes the grid, the problem, and
// the solution vector. It then sets up the time loop and
// the VTK output module. The main loop advances the solution
// in time and writes the solution to VTK files.
//
// The assembler assembles the linear system in residual form
// and the linear solver solves the linear system. Here, we
// use an algebraic multigrid preconditioned conjugate gradient
// (AMG-CG) solver.
//
///////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using TypeTag = Properties::TTag::DiffusionTest;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    GridManager<Grid> gridManager;
    gridManager.init();
    const auto gridData = gridManager.getGridData();

    const auto leafGridView = gridManager.grid().leafGridView();
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    auto problem = std::make_shared<Problem>(gridGeometry);
    SolutionVector sol(gridGeometry->numDofs()); sol = 0.0;
    SolutionVector dataSol = sol; // for output

    // read diffusion tensor data and concentration boundary data
    problem->setGridData(*gridData);

    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(
        Chrono::toSeconds(getParam("TimeLoop.TStart", "0")),
        Chrono::toSeconds(getParam("TimeLoop.InitialTimeStepSize")),
        Chrono::toSeconds(getParam("TimeLoop.TEnd"))
    );

    // Concentration averages for different regions
    // assume markers start consecutive from 1, ... and 0 is the quantity for the whole mesh
    constexpr std::size_t numRegions = 4;
    std::array<std::string, numRegions> regions{"total", "corticalgraymatter", "whitematter", "subcorticalgraymatter"};

    // data arrays for outout data
    std::array<std::vector<Dune::FieldVector<double, 2>>, numRegions>
        concentrationData, amountData, concentrationDataMRI, amountDataMRI;
    if (gridGeometry->gridView().comm().rank() == 0)
    {
        for (int i = 0; i < numRegions; ++i)
        {
            concentrationData[i].reserve(72);
            amountData[i].reserve(72);
            concentrationDataMRI[i].reserve(5);
            amountDataMRI[i].reserve(5);
        }
    }

    // function to compute output data points
    const auto computeRegionalAverages = [&](bool outputVolumeSummary = false)
    {
        const auto& gridView = gridGeometry->gridView();
        const bool isCheckPoint = timeLoop->isCheckPoint() || timeLoop->timeStepIndex() == 0;
        if (isCheckPoint)
        {
            // read the grid data into a solution vector
            const int t = timeLoop->time();
            for (const auto& vertex : vertices(gridView))
            {
                const auto vIdx = gridGeometry->dofMapper().index(vertex);
                dataSol[vIdx] = gridData->getParameter(vertex, Fmt::format("concentration-{:d}", t));
            }
        }

        std::array<Scalar, numRegions> amount = {}, amountMRI = {}, volume = {};
        for (const auto& element : elements(gridView, Dune::Partitions::interior))
        {
            const auto elemSol = elementSolution(element, sol, *gridGeometry);
            const auto geometry = element.geometry();
            const auto center = geometry.center();
            const auto vol = geometry.volume();
            const auto c = evalSolution(element, geometry, *gridGeometry, elemSol, center);
            amount[0] += c*vol;
            volume[0] += vol;

            Scalar cMRI = 0.0;
            if (isCheckPoint)
            {
                const auto elemSol = elementSolution(element, dataSol, *gridGeometry);
                cMRI = evalSolution(element, geometry, *gridGeometry, elemSol, center);
                amountMRI[0] += cMRI*vol;
            }

            const auto marker = gridData->getParameter(element, "subdomains");
            if (!(marker == 1 || marker == 2 || marker == 3))
                DUNE_THROW(Dune::IOError, "Unknown subdomain marker " << marker);

            amount[marker] += c*vol;
            volume[marker] += vol;
            amountMRI[marker] += cMRI*vol;
        }

        for (int i = 0; i < numRegions; ++i)
        {
            amount[i] = gridView.comm().sum(amount[i]);
            volume[i] = gridView.comm().sum(volume[i]);
        }

        if (isCheckPoint)
            for (int i = 0; i < numRegions; ++i)
                amountMRI[i] = gridView.comm().sum(amountMRI[i]);

        if (gridView.comm().rank() == 0)
        {
            const Scalar l_per_mm3 = 1.0e6; // mesh is in mm and concentration in mol/m^3
            const Scalar t = timeLoop->time();

            for (int i = 0; i < numRegions; ++i)
            {
                concentrationData[i].push_back({t, amount[i]/volume[i]}); // in mol/m^3 = mmol/L
                amountData[i].push_back({t, amount[i]/l_per_mm3}); // in mmol
            }

            if (isCheckPoint)
            {
                for (int i = 0; i < numRegions; ++i)
                {
                    concentrationDataMRI[i].push_back({t, amountMRI[i]/volume[i]}); // in mol/m^3 = mmol/L
                    amountDataMRI[i].push_back({t, amountMRI[i]/l_per_mm3}); // in mmol
                }
            }

            if (outputVolumeSummary)
            {
                for (int i = 0; i < numRegions; ++i)
                {
                    std::cout << Fmt::format("Region {}\n", regions[i]);
                    std::cout << Fmt::format("Volume {}: {} mm^3\n", regions[i], volume[i]);
                    std::cout << Fmt::format("Amount {}: {} mmol\n", regions[i], amount[i]/l_per_mm3);
                    std::cout << Fmt::format("Concentration {}: {} mmol/L\n", regions[i], amount[i]/volume[i]);
                }
            }
        }
    };

    computeRegionalAverages(/*outputVolumeSummary=*/true);

    // VTK output
    IO::OutputModule vtkWriter(
        IO::Format::pvd_with(
            IO::Format::vtu.with({
                .encoder = IO::Encoding::raw,
                .compressor = IO::Compression::none,
                .data_format = IO::VTK::DataFormat::appended
            })
        ),
        *gridVariables, sol, problem->name()
    );

    vtkWriter.addVolumeVariable([](const auto& vv){ return vv.priVar(0); }, "c");
    std::vector<double> boundaryConcentration(gridGeometry->numDofs(), 0.0);
    problem->updateBoundaryDataForOutput(boundaryConcentration, timeLoop->time());
    vtkWriter.addField(boundaryConcentration, "c_b");
    vtkWriter.addPointField([&](const auto& v){ return dataSol[gridGeometry->dofMapper().index(v)][0]; }, "c_data");
    vtkWriter.write(0.0);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = AMGCGIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using Solver = Dumux::LinearPDESolver<Assembler, LinearSolver>;

    auto oldSol = sol; // copy the vector to store state of previous time step
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, oldSol);
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
    Solver solver(assembler, linearSolver);

    const bool vtkAtCheckPoints = getParam<bool>("VTK.WriteOnlyAtCheckPoints", true);
    const auto& cp = problem->checkPoints();
    timeLoop->setCheckPoint(cp.begin()+1, cp.end());
    timeLoop->start(); do
    {
        // set the current time in the problem
        problem->setTime(timeLoop->time() + timeLoop->timeStepSize());

        // assemble & solve
        solver.solve(sol);

        // make the new solution the old solution
        oldSol = sol;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write VTK output (writes out the concentration field)
        problem->updateBoundaryDataForOutput(boundaryConcentration, timeLoop->time());
        computeRegionalAverages();

        if (vtkAtCheckPoints && timeLoop->isCheckPoint())
            vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

    } while (!timeLoop->finished());

    timeLoop->finalize(gridGeometry->gridView().comm());

    if (gridGeometry->gridView().comm().rank() == 0)
    {
        for (int i = 0; i < numRegions; ++i)
        {
            writeContainerToFile(concentrationData[i], Fmt::format("{}.txt", regions[i]));
            writeContainerToFile(amountData[i], Fmt::format("{}_amount.txt", regions[i]));
            writeContainerToFile(concentrationDataMRI[i], Fmt::format("{}_data.txt", regions[i]));
            writeContainerToFile(amountDataMRI[i], Fmt::format("{}_amount_data.txt", regions[i]));
        }
    }

    return 0;
}
