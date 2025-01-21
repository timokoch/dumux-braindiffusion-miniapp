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

namespace Dumux::Properties::TTag {

// compile-time configuration of simulation
struct DiffusionTest
{
    using InheritsFrom = std::tuple<DiffusionModel, BoxModel>;

    using Scalar = double;
    using Grid = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

} // end namespace Dumux::Properties::TTag

double concentrationBoundaryDataFunction_(const double time, const Dune::BlockVector<double>& p)
{
    // check if any parameters are not finite
    for (int i = 0; i < p.size(); ++i)
        if (!std::isfinite(p[i]))
            DUNE_THROW(Dune::Exception, "Invalid parameter values: " << p);

    const auto days = time/86400.0;
    return p[0]*days*std::exp(-p[1]*days);
    // return p[0]*(std::exp(-p[1]*days) - std::exp(-p[2]*days));
}

//////////////////////////////////////////////////////////////
//
// Compute curve fits for the concentration boundary data.
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

    GridManager<Grid> gridManager;
    gridManager.init();
    const auto gridData = gridManager.getGridData();

    const auto leafGridView = gridManager.grid().leafGridView();
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    const auto& gridView = gridGeometry->gridView();

    Dune::BlockVector<Scalar> times = {0, 15716, 86001, 172280, 251641};
    const auto plotTimes = linspace<double>(0, 251641, 100);
    auto plotCValues = plotTimes;

    writeContainerToFile(times, "times.dat");
    writeContainerToFile(plotTimes, "plottimes.dat");

    Dune::BlockVector<Scalar> boundaryData(times.size());
    std::vector<Dune::BlockVector<Scalar>> concentrationBoundaryDataParams(gridGeometry->numDofs());
    std::size_t converged = 0;
    std::size_t numBoundaryVertices = 0;
    for (const auto& vertex : vertices(gridView))
    {
        const auto vIdx = gridGeometry->dofMapper().index(vertex);
        if (!gridGeometry->dofOnBoundary(vIdx))
            continue;

        ++numBoundaryVertices;

        for (int i = 0; i < times.size(); ++i)
        {
            boundaryData[i]
                = gridData->getParameter(vertex, Fmt::format("boundary-concentration-{:d}", int(times[i])));

            if (!std::isfinite(boundaryData[i]))
                DUNE_THROW(Dune::Exception, "Invalid boundary data at vertex " << vIdx << " at time " << times[i]);
        }

        Dune::BlockVector<double> p({0.1, 1.0});
        // Dune::BlockVector<double> p({0.1, 1.0, 2.0});
        Dune::BlockVector<double> res(times.size());
        const auto residual = [&](const auto& params)
        {
            for (size_t i = 0; i < times.size(); ++i)
            {
                res[i] = concentrationBoundaryDataFunction_(times[i], params) - boundaryData[i];

                // enforce some bounds on the parameters
                if (params[0] < 0.0 || params[0] > 10.0)
                    res[i] = params[0]*1e4;

                if (params[1] < 0.0 || params[1] > 5.0)
                    res[i] = params[1]*1e4;

                if (params.size() > 2)
                {
                    if (params[2] < 0.0 || params[2] > 10.0)
                        res[i] = params[2]*1e4;

                    if (std::abs(params[1] - params[2]) < 1e-2)
                        res[i] = std::max(1e4, 1.0/(params[1] - params[2])*1e4);

                    if (params[1] > params[2])
                        res[i] = params[1]*params[2]*1e4;
                }

                if (!std::isfinite(res[i]))
                    DUNE_THROW(Dune::Exception, "Invalid residual value at vertex " << vIdx << " at time "
                                << times[i] << ": " << res[i] << " for parameters " << params);
            }
            return res;
        };

        {

            std::ostringstream local;
            auto cout_buff = std::cout.rdbuf();
            std::cout.rdbuf(local.rdbuf());

            const auto c = Optimization::makeNonlinearLeastSquaresSolver(residual, p, times.size())->apply(p);

            std::cout.rdbuf(cout_buff);

            if (c) ++converged;
        }

        for (int i = 0; i < plotTimes.size(); ++i)
            plotCValues[i] = concentrationBoundaryDataFunction_(plotTimes[i], p);

        writeContainerToFile(plotCValues, Fmt::format("curvefit-{:d}.dat", vIdx));
        writeContainerToFile(boundaryData, Fmt::format("data-{:d}.dat", vIdx));
    }

    std::cout << "Curve fit converged for " << converged
              << " / " << numBoundaryVertices << std::endl;

    return 0;
}
