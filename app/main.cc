// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
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

#include <dumux/io/vtkoutputmodule.hh>
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
        return concentrationBoundaryDataFunction(time_, p);
    }

    Scalar concentrationBoundaryDataFunction(const Scalar time, const Dune::BlockVector<Scalar>& p) const
    {
        // check if any parameters are not finite
        for (int i = 0; i < p.size(); ++i)
            if (!std::isfinite(p[i]))
                DUNE_THROW(Dune::Exception, "Invalid parameter values: " << p);

        if (useCurveFit_)
            return p[0]*(std::exp(std::clamp(p[1], 0.01, 10.0)*time/86400.0) - std::exp(-std::clamp(p[2], 0.01, 10.0)*time/86400.0));
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
                Dune::BlockVector<double> p({0.2, 2.0, 0.2});
                Dune::BlockVector<double> res(times_.size());
                const auto residual = [&](const auto& params)
                {
                    for (size_t i = 0; i < times_.size(); ++i)
                    {
                        res[i] = this->concentrationBoundaryDataFunction(times_[i], params) - boundaryData[i];
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
            std::cout << "Biexponential curve fit converged for " << converged
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
    std::vector<Dune::FieldVector<double, 2>> totalConcentrationData, grayMatterConcentrationData, whiteMatterConcentrationData;
    std::vector<Dune::FieldVector<double, 2>> totalAmountData, grayMatterAmountData, whiteMatterAmountData;
    std::vector<Dune::FieldVector<double, 2>> totalConcentrationDataData, grayMatterConcentrationDataData, whiteMatterConcentrationDataData;
    std::vector<Dune::FieldVector<double, 2>> totalAmountDataData, grayMatterAmountDataData, whiteMatterAmountDataData;
    if (gridGeometry->gridView().comm().rank() == 0)
    {
        totalConcentrationData.reserve(72);
        grayMatterConcentrationData.reserve(72);
        whiteMatterConcentrationData.reserve(72);
        totalAmountData.reserve(72);
        grayMatterAmountData.reserve(72);
        whiteMatterAmountData.reserve(72);
        totalConcentrationDataData.reserve(5);
        grayMatterConcentrationDataData.reserve(5);
        whiteMatterConcentrationDataData.reserve(5);
        totalAmountDataData.reserve(5);
        grayMatterAmountDataData.reserve(5);
        whiteMatterAmountDataData.reserve(5);
    }

    const auto computeRegionalAverages = [&]()
    {
        const auto& gridView = gridGeometry->gridView();

        Scalar totalAmount = 0.0, grayMatterAmount = 0.0, whiteMatterAmount = 0.0;
        Scalar totalAmountInData = 0.0, grayMatterAmountInData = 0.0, whiteMatterAmountInData = 0.0;
        Scalar totalVolume = 0.0, grayMatterVolume = 0.0, whiteMatterVolume = 0.0;

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

        for (const auto& element : elements(gridView, Dune::Partitions::interior))
        {
            const auto elemSol = elementSolution(element, sol, *gridGeometry);
            const auto geometry = element.geometry();
            const auto center = geometry.center();
            const auto volume = geometry.volume();
            const auto c = evalSolution(element, geometry, *gridGeometry, elemSol, center);
            totalAmount += c*volume;
            totalVolume += volume;

            Scalar cData = 0.0;
            if (isCheckPoint)
            {
                const auto elemSol = elementSolution(element, dataSol, *gridGeometry);
                cData = evalSolution(element, geometry, *gridGeometry, elemSol, center);
                totalAmountInData += cData*volume;
            }

            const auto marker = gridData->getParameter(element, "subdomains");
            if (marker == 1)
            {
                grayMatterAmount += c*volume;
                grayMatterVolume += volume;
                grayMatterAmountInData += cData*volume;
            }
            else if (marker == 2)
            {
                whiteMatterAmount += c*volume;
                whiteMatterVolume += volume;
                whiteMatterAmountInData += cData*volume;
            }
        }

        totalAmount = gridView.comm().sum(totalAmount);
        grayMatterAmount = gridView.comm().sum(grayMatterAmount);
        whiteMatterAmount = gridView.comm().sum(whiteMatterAmount);
        grayMatterVolume = gridView.comm().sum(grayMatterVolume);
        whiteMatterVolume = gridView.comm().sum(whiteMatterVolume);
        totalVolume = gridView.comm().sum(totalVolume);

        if (isCheckPoint)
        {
            totalAmountInData = gridView.comm().sum(totalAmountInData);
            grayMatterAmountInData = gridView.comm().sum(grayMatterAmountInData);
            whiteMatterAmountInData = gridView.comm().sum(whiteMatterAmountInData);
        }

        if (gridView.comm().rank() == 0)
        {
            const auto t = timeLoop->time();
            totalConcentrationData.push_back({t, totalAmount/totalVolume});
            grayMatterConcentrationData.push_back({t, grayMatterAmount/grayMatterVolume});
            whiteMatterConcentrationData.push_back({t, whiteMatterAmount/whiteMatterVolume});
            totalAmountData.push_back({t, totalAmount});
            grayMatterAmountData.push_back({t, grayMatterAmount});
            whiteMatterAmountData.push_back({t, whiteMatterAmount});

            if (isCheckPoint)
            {
                totalConcentrationDataData.push_back({t, totalAmountInData/totalVolume});
                grayMatterConcentrationDataData.push_back({t, grayMatterAmountInData/grayMatterVolume});
                whiteMatterConcentrationDataData.push_back({t, whiteMatterAmountInData/whiteMatterVolume});
                totalAmountDataData.push_back({t, totalAmountInData});
                grayMatterAmountDataData.push_back({t, grayMatterAmountInData});
                whiteMatterAmountDataData.push_back({t, whiteMatterAmountInData});
            }
        }
    };

    computeRegionalAverages();

    // VTK output
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, sol, problem->name());
    vtkWriter.addVolumeVariable([](const auto& vv){ return vv.priVar(0); }, "c");
    std::vector<double> boundaryConcentration(gridGeometry->numDofs(), 0.0);
    problem->updateBoundaryDataForOutput(boundaryConcentration, timeLoop->time());
    vtkWriter.addField(boundaryConcentration, "c_b");
    vtkWriter.addField(dataSol, "c_data");
    vtkWriter.write(0.0);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = AMGCGIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using Solver = Dumux::LinearPDESolver<Assembler, LinearSolver>;

    auto oldSol = sol; // copy the vector to store state of previous time step
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, oldSol);
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
    Solver solver(assembler, linearSolver);

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
        vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

    } while (!timeLoop->finished());

    timeLoop->finalize(gridGeometry->gridView().comm());

    if (gridGeometry->gridView().comm().rank() == 0)
    {
        writeContainerToFile(totalConcentrationData, "total.txt");
        writeContainerToFile(grayMatterConcentrationData, "graymatter.txt");
        writeContainerToFile(whiteMatterConcentrationData, "whitematter.txt");
        writeContainerToFile(totalAmountData, "total_amount.txt");
        writeContainerToFile(grayMatterAmountData, "graymatter_amount.txt");
        writeContainerToFile(whiteMatterAmountData, "whitematter_amount.txt");

        writeContainerToFile(totalConcentrationDataData, "total_data.txt");
        writeContainerToFile(grayMatterConcentrationDataData, "graymatter_data.txt");
        writeContainerToFile(whiteMatterConcentrationDataData, "whitematter_data.txt");
        writeContainerToFile(totalAmountDataData, "total_amount_data.txt");
        writeContainerToFile(grayMatterAmountDataData, "graymatter_amount_data.txt");
        writeContainerToFile(whiteMatterAmountDataData, "whitematter_amount_data.txt");
    }

    return 0;
}
