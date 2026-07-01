using System;
using System.IO;
using System.Text;
using Archivory;

namespace Clustering.Graphical.SPC.Runtime.Core.Sampler
{
    /// <summary>
    /// Binary serializer for <see cref="Accumulator"/>, owning the SPCX checkpoint format. Access via
    /// the static <see cref="Instance"/> singleton:
    /// <code>
    /// AccumulatorSerializer.Instance.WriteToFile(accumulator, path);
    /// </code>
    /// </summary>
    /// <remarks>
    /// Little-endian, 4-byte magic + version. Writes the scalar moments + resume state, then — when the
    /// run tracked them — the per-edge arrays inline behind a presence flag. The old paired <c>.spce</c>
    /// sidecar and its key-matching are retired: the per-edge section is now positionally nested in the
    /// same file (one accumulator, one artifact). Composes on Archivory's
    /// <see cref="BinarySerializerBase{T}"/> for atomic file I/O — this class supplies only the
    /// SPC-specific stream schema.
    /// </remarks>
    public sealed class AccumulatorSerializer : BinarySerializerBase<Accumulator>
    {
        public static AccumulatorSerializer Instance { get; } = new();
        public override string DefaultFileExtension => ".spcx";

        private const uint Magic = 0x58435053; // "SPCX"
        private const int Version = 7;         // v7: add co-membership count array behind its own presence flag

        public override void WriteTo(Accumulator value, Stream stream)
        {
            if (value is null)
                throw new ArgumentNullException(nameof(value));

            using var w = new BinaryWriter(stream, Encoding.UTF8, leaveOpen: true);
            w.Write(Magic);
            w.Write(Version);
            w.Write(value.Temperature);
            w.Write(value.Q);
            w.Write(value.ReplicaIndex);
            w.Write(value.DrawCount);

            w.Write(value.Spins.Length);
            foreach (int s in value.Spins) w.Write(s);

            w.Write(value.ClusterSizeHistogram.Length);
            foreach (int h in value.ClusterSizeHistogram) w.Write(h);

            w.Write(value.RngState0);
            w.Write(value.RngState1);
            w.Write(value.RngState2);
            w.Write(value.RngState3);

            w.Write(value.RunningSumSqClusterSizes);
            w.Write(value.RunningSumSqClusterSizesExcl);
            w.Write(value.RunningSumEnergy);
            w.Write(value.RunningSumEnergySq);
            w.Write(value.RunningSumMag);
            w.Write(value.RunningSumMagSq);

            // Per-edge section (v7: independently-gated flags — v6 used one shared flag which broke
            // when Affinities=true + Alignments=false).
            bool hasBond = value.BondFormedCount is not null;
            w.Write(hasBond);
            if (hasBond)
            {
                w.Write(value.BondFormedCount!.Length);
                foreach (int x in value.BondFormedCount) w.Write(x);
            }

            bool hasSpin = value.SpinAgreementCount is not null;
            w.Write(hasSpin);
            if (hasSpin)
            {
                w.Write(value.SpinAgreementCount!.Length);
                foreach (int x in value.SpinAgreementCount) w.Write(x);
            }

            // Per-node section (v6): the two landscape sums are independently gated,
            // so each carries its own presence flag + length (mirrors the per-edge nesting).
            bool hasNodeSize = value.SumClusterSizePerNode is not null;
            w.Write(hasNodeSize);
            if (hasNodeSize)
            {
                w.Write(value.SumClusterSizePerNode!.Length);
                foreach (double x in value.SumClusterSizePerNode) w.Write(x);
            }

            bool hasNodeGiant = value.SumInGiantClusterPerNode is not null;
            w.Write(hasNodeGiant);
            if (hasNodeGiant)
            {
                w.Write(value.SumInGiantClusterPerNode!.Length);
                foreach (double x in value.SumInGiantClusterPerNode) w.Write(x);
            }

            // Co-membership section (v7): presence flag + edge-length array when tracked.
            bool hasCoMembership = value.CoMembershipCount is not null;
            w.Write(hasCoMembership);
            if (hasCoMembership)
            {
                w.Write(value.CoMembershipCount!.Length);
                foreach (int x in value.CoMembershipCount) w.Write(x);
            }

            w.Write(value.TimestampUtc.Ticks);
        }

        public override Accumulator ReadFrom(Stream stream)
        {
            using var r = new BinaryReader(stream, Encoding.UTF8, leaveOpen: true);

            uint magic = r.ReadUInt32();
            if (magic != Magic)
                throw new InvalidDataException(
                    $"Not an SPCX checkpoint (magic 0x{magic:X8}, expected 0x{Magic:X8}).");
            int version = r.ReadInt32();
            if (version != Version)
                throw new InvalidDataException(
                    $"Unsupported SPCX version {version} (this reader handles version {Version}).");

            double temperature  = r.ReadDouble();
            int    q            = r.ReadInt32();
            int    replicaIndex = r.ReadInt32();
            int    drawCount    = r.ReadInt32();

            int spinCount = r.ReadInt32();
            var spins = new int[spinCount];
            for (int i = 0; i < spinCount; i++) spins[i] = r.ReadInt32();

            int histBins = r.ReadInt32();
            var hist = new int[histBins];
            for (int i = 0; i < histBins; i++) hist[i] = r.ReadInt32();

            ulong s0 = r.ReadUInt64();
            ulong s1 = r.ReadUInt64();
            ulong s2 = r.ReadUInt64();
            ulong s3 = r.ReadUInt64();

            double sumSqCluster     = r.ReadDouble();
            double sumSqClusterExcl = r.ReadDouble();
            double sumEnergy        = r.ReadDouble();
            double sumEnergySq      = r.ReadDouble();
            double sumMag           = r.ReadDouble();
            double sumMagSq         = r.ReadDouble();

            int[]? bond = null;
            bool hasBond = r.ReadBoolean();
            if (hasBond)
            {
                int edgeCount = r.ReadInt32();
                bond = new int[edgeCount];
                for (int i = 0; i < edgeCount; i++) bond[i] = r.ReadInt32();
            }

            int[]? spin = null;
            bool hasSpin = r.ReadBoolean();
            if (hasSpin)
            {
                int edgeCount = r.ReadInt32();
                spin = new int[edgeCount];
                for (int i = 0; i < edgeCount; i++) spin[i] = r.ReadInt32();
            }

            // Per-node section (v6): two independently-gated landscape sums.
            double[]? nodeSize = null;
            bool hasNodeSize = r.ReadBoolean();
            if (hasNodeSize)
            {
                int len = r.ReadInt32();
                nodeSize = new double[len];
                for (int i = 0; i < len; i++) nodeSize[i] = r.ReadDouble();
            }

            double[]? nodeGiant = null;
            bool hasNodeGiant = r.ReadBoolean();
            if (hasNodeGiant)
            {
                int len = r.ReadInt32();
                nodeGiant = new double[len];
                for (int i = 0; i < len; i++) nodeGiant[i] = r.ReadDouble();
            }

            // Co-membership section (v7).
            int[]? coMembership = null;
            bool hasCoMembership = r.ReadBoolean();
            if (hasCoMembership)
            {
                int len = r.ReadInt32();
                coMembership = new int[len];
                for (int i = 0; i < len; i++) coMembership[i] = r.ReadInt32();
            }

            long ticks = r.ReadInt64();

            return new Accumulator
            {
                Temperature              = temperature,
                Q                        = q,
                ReplicaIndex             = replicaIndex,
                DrawCount                = drawCount,
                Spins                    = spins,
                ClusterSizeHistogram     = hist,
                RngState0                = s0,
                RngState1                = s1,
                RngState2                = s2,
                RngState3                = s3,
                RunningSumSqClusterSizes     = sumSqCluster,
                RunningSumSqClusterSizesExcl = sumSqClusterExcl,
                RunningSumEnergy             = sumEnergy,
                RunningSumEnergySq           = sumEnergySq,
                RunningSumMag                = sumMag,
                RunningSumMagSq              = sumMagSq,
                BondFormedCount    = bond,
                SpinAgreementCount = spin,
                CoMembershipCount  = coMembership,
                SumClusterSizePerNode    = nodeSize,
                SumInGiantClusterPerNode = nodeGiant,
                TimestampUtc       = new DateTime(ticks, DateTimeKind.Utc),
            };
        }
    }
}
