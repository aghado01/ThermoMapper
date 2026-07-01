#nullable enable
using System;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// A persistent topological cell used primarily as the currency for Zigzag Persistent Homology.
/// Represents an immutable n-dimensional entity with a permanent global ID and a specific Z/2 boundary.
/// </summary>
public interface IAbstractCell
{
    /// <summary>
    /// Permanent global cell identifier, assigned once and never reused.
    /// </summary>
    int GlobalCellId { get; }

    /// <summary>
    /// Dimension of the cell (0 = vertex, 1 = edge, 2 = face, etc).
    /// </summary>
    int Dimension { get; }

    /// <summary>
    /// Z/2 boundary as an array of indices into earlier (p-1)-cells.
    /// Note: Indices refer to cells already present in the filtration immediately before this cell's addition.
    /// </summary>
    int[] Boundary { get; }
}

/// <summary>
/// Helper for deriving face indices against a simplicial complex.
/// </summary>
public static class SimplicialCell
{
    /// <summary>
    /// Derives face indices from the vertex set against the cells already present.
    /// Used to generate boundary arrays for <see cref="IAbstractCell"/> implementors.
    /// </summary>
    public static int[] GetBoundaryIndices(ReadOnlySpan<int> vertices, IFiltration currentComplex)
    {
        // For a simplex with vertices v0, v1, ..., vn, the boundary faces are obtained
        // by omitting one vertex at a time. This method maps those sub-simplices to their
        // indices in the current complex.
        
        if (vertices.Length <= 1)
        {
            return Array.Empty<int>(); // Vertices have empty boundaries
        }

        int[] boundaryIndices = new int[vertices.Length];
        Span<int> faceVertices = stackalloc int[vertices.Length - 1];

        for (int i = 0; i < vertices.Length; i++)
        {
            // Omit the i-th vertex
            int dst = 0;
            for (int j = 0; j < vertices.Length; j++)
            {
                if (i != j)
                {
                    faceVertices[dst++] = vertices[j];
                }
            }

            // Find the face in the current complex
            int faceIndex = -1;
            for (int k = 0; k < currentComplex.Count; k++)
            {
                if (currentComplex.GetDimension(k) == vertices.Length - 2)
                {
                    var cVertices = currentComplex.GetVertices(k);
                    bool match = true;
                    // Assumes vertices are sorted/comparable in the same way
                    for (int v = 0; v < faceVertices.Length; v++)
                    {
                        if (cVertices[v] != faceVertices[v])
                        {
                            match = false;
                            break;
                        }
                    }

                    if (match)
                    {
                        faceIndex = k;
                        break;
                    }
                }
            }

            if (faceIndex == -1)
            {
                throw new InvalidOperationException("Face not found in the current complex. Ensure the complex is built correctly.");
            }

            boundaryIndices[i] = faceIndex;
        }

        return boundaryIndices;
    }
}
