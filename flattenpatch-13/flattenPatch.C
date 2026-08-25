/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    flattenPatch

Description
    Flattens the patch that is specified in the dictionary to the set direction
    and coordinate position

Author
    Tobias Holzmann

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "SubField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    Info<< "Read flattenPatchDict\n" << endl;

    IOdictionary flattenPatchDict
    (
        IOobject
        (
            "flattenPatchDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    //- Get the patch names
    const dictionary& patchesToBeManipulated =
        flattenPatchDict.subDict("patches");

    const tokenList patchName = patchesToBeManipulated.tokens();
    wordList patchNames;

    //- Get normal direction
    const word vecDir = flattenPatchDict.lookup("direction");

    //- Initialize with x direction
    direction normalDir = vector::X;

    if (vecDir == "x" || vecDir == "X")
    {
        normalDir = vector::X;
    }
    else if (vecDir == "y" || vecDir == "Y")
    {
        normalDir = vector::Y;
    }
    else if (vecDir == "z" || vecDir == "Z")
    {
        normalDir = vector::Z;
    }
    else
    {
        FatalErrorInFunction
            << " Direction of the vector is not correct. Available options are"
            << "\n x\n y\n z" << nl << endl
            << abort(FatalError);
    }

    //- Get Cmp value
    const scalar cmpVal = readScalar(flattenPatchDict.lookup("setCmpTo"));

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(polyMesh::meshSubDir, "points"),
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    //- Algorithm
    //  Find the patch which should be modified and move the points
    //  corresponding to the normal vector we set
    //  If not found, do nothing

    const polyBoundaryMesh& boundaryPatches = mesh.boundaryMesh();
    wordList patchModified;


    forAll(patchName, i)
    {
        if (patchName[i].isWord())
        {
            bool found = false;

            const polyPatch* boundaryPatch = NULL;

            forAll(boundaryPatches, patchi)
            {
                const token tmp = boundaryPatches[patchi].name();

                if (tmp == patchName[i])
                {
                    found = true;
                    boundaryPatch = &boundaryPatches[patchi];
                    break;
                }
            }

            if (!found)
            {
                FatalErrorInFunction
                    << " Patch '" << patchName[i] << "' not found.\n"
                    << " Valid patches are:\n " << boundaryPatches.names()
                    << abort(FatalError);
            }
            else
            {
                patchModified.append(boundaryPatch->name());

                //- Patch normas
                const vectorField::subField pN = boundaryPatch->faceAreas();

                //- Averaged pN
                const vector avePN = gAverage(pN);
                const vector magAvePN = avePN / mag(avePN);

                Info<< "Average patch normal of patch '"
                    << boundaryPatch->name() << "' is " << avePN << nl;

                Info<< "Normalized average patch normal of patch '"
                    << boundaryPatch->name() << "' is " << magAvePN << nl;

                //- Get point IDs of the patch
                const labelList& boundaryPointID = boundaryPatch->meshPoints();

                Info<< "Modifing the points on the patch\n" << nl;

                forAll(boundaryPointID, celli)
                {
                    const label ID = boundaryPointID[celli];

                    //- Set the corresponding point value
                    points[ID][normalDir] = cmpVal;
                }
            }
        }
    }

    if (patchModified.empty())
    {
        Pout<< "No patch was modified\n" << endl;
    }
    else
    {
        forAll(patchModified, i)
        {
            Info<< "Patch '" << patchModified[i] << "' was modified\n" << endl;
        }
    }

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< "Writing points into directory " << points.path() << nl << endl;
    points.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
