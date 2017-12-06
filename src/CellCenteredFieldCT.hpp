#pragma once
#include "Array.hpp"
#include "Mara.hpp"




/**
A class to encapsulate a limited set of algorithms for manipulating cell-
centered magnetic field operationrs. These operations do not account for any
mesh geometry; the array data provided is assumed to be logically cartesian
and uniformly spaced with lattice spacing equal to meshSpacing on each axis,
which by default is set to 1.0.
*/
class CellCenteredFieldCT
{
public:
    using Array = Cow::Array;

    /** Constructor */
    CellCenteredFieldCT();

    /** Use the given mesh spacing for finite difference operations. */
    void setMeshSpacing (double uniformMeshSpacing) { meshSpacing = uniformMeshSpacing; }

    /**
    Change the input array of face-centered (Godunov) fluxes to yield
    divergenceless magnetic fluxes. The array may contain fluxes of other
    quantities on axis 3, but the fluxes in the three cardinal directions must
    be contiguous, and begin at magneticIndex.
    */
    void correctGodunovFluxes (Array& F, int magneticIndex) const;

    /**
    Same as correctGodunovFluxes, except a new array is returned instead of
    doing the correction in-place.
    */
    Array generateGodunovFluxes (const Array& F, int magneticIndex) const;

    /**
    Given an array of vector potential (or electric field) with shape (nx, ny,
    nz, 3, 3), where axis 3 labels the vector component, and axis 4 labels the
    face, return ideal MHD fluxes, i.e. Fi(Bj) = -Ak where i,j,k are a
    positive permutation of 1,2,3. The array returned has the same shape as B.
    */
    Array vectorPotentialToFluxes (Array A) const;

    /**
    Return the field-CT, vertex-situated div.B, where B is the cell-centered
    magnetic field. If location=MeshLocation::vert then the 8 cells sharing a
    vertex are used to compute the divergence, and the returned array has the
    shape of B, increased by 1 on each axis. If the mesh location is cell,
    then the returned array has the same shape as B, and the divergence at
    each cell's 6 vertices are averaged to estimate the divergence at the cell
    center.
    */
    Array monopole (Array B, MeshLocation location) const;

    /**
    Return an estimate of curl.B, where B is the cell-centered magnetic field.
    Presently, only location=MeshLocation::cell is supported; a centered
    second-order derivative operator is utilized to compute the gradient, and
    the current on the outer-most layer of cells (1 deep) is left zero, but
    the returned array has the same shape as B.
    */
    Array current (Array B, MeshLocation location) const;

private:
    double meshSpacing;
};
