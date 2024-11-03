// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_SIMPLE_DIFFUSION_MODEL_HH
#define DUMUX_SIMPLE_DIFFUSION_MODEL_HH

//////////////////////////////////////////////////////////////
//
// A simple diffusion model:
//
//     ∂c/∂t - ∇·(D∇c) = q,
//
// where c is total concentration, D is the diffusion tensor,
// and q is a source term. The source term is implemented
// in the `source` method of the problem class.
//
// The equations are discreized with a control-volume finite
// element method with piece-wise linear shape functions on
// elements (P1).
//
///////////////////////////////////////////////////////////////

#include <dune/common/fvector.hh>
#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/volumevariables.hh>
#include <dumux/discretization/method.hh>

namespace Dumux::Properties::TTag {
struct DiffusionModel {};
} // end namespace Dumux::Properties:TTag

namespace Dumux {

template<class TypeTag>
class DiffusionModelLocalResidual
: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = typename GridVariables::GridVolumeVariables::VolumeVariables;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    using ParentType::ParentType;

    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage;
        storage[Indices::massBalanceEqIdx]
            = volVars.priVar(Indices::concentrationIdx);
        return storage;
    }

    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        static_assert(DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>,
            "This local residual is hard-coded to control-volume finite element schemes");

        // Compute ∇c at the integration point of this sub control volume face.
        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, dimWorld> gradConcentration(0.0);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            gradConcentration.axpy(
                volVars.priVar(Indices::concentrationIdx),
                fluxVarCache.gradN(scv.indexInElement())
            );
        }

        NumEqVector flux;

        // Compute the flux with `vtmv` (vector transposed times matrix times vector) or -n^T D ∇c A.
        flux[Indices::massBalanceEqIdx] = -1.0*vtmv(
            scvf.unitOuterNormal(), problem.diffusionCoefficient(element, fvGeometry.scv(scvf.insideScvIdx())), gradConcentration
        )*scvf.area();

        return flux;
    }
};

} // end namespace Dumux

namespace Dumux::Properties {

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::DiffusionModel>
{ using type = DiffusionModelLocalResidual<TypeTag>; };

template<class TypeTag>
struct Scalar<TypeTag, TTag::DiffusionModel>
{ using type = double; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::DiffusionModel>
{
    struct type
    {
        struct Indices
        {
            static constexpr int concentrationIdx = 0;
            static constexpr int massBalanceEqIdx = 0;
        };

        static constexpr int numEq() { return 1; }
    };
};

template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::DiffusionModel>
{
    using type = Dune::FieldVector<
        GetPropType<TypeTag, Properties::Scalar>,
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()
    >;
};

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::DiffusionModel>
{
    struct Traits
    {
        using PrimaryVariables
            = GetPropType<TypeTag, Properties::PrimaryVariables>;
    };
    using type = BasicVolumeVariables<Traits>;
};

} // end namespace Dumux::Properties

#endif
