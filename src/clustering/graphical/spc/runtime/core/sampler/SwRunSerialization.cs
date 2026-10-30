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
    /// Little-endian, 4-byte magic + version. Writes the scalar moments + resume state and the per-draw
    /// record, then — when the run tracked them — the per-edge arrays inline behind a presence flag. The old paired <c>.spce</c>
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
        private const int Version = 11;        // v11: the exact cluster-size histogram (v10: the burn-in trace; v9: the per-draw record)

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
            foreach (long h in value.ClusterSizeHistogram) w.Write(h);

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

            // Per-draw record (v9): always present.
            WriteRecord(w, "Per-draw record",
                value.EnergyPerDraw, value.MagPerDraw, value.SumSqClusterSizesPerDraw, value.MaxClusterSizePerDraw);

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

            // Burn-in trace (v10): presence flag, then the per-draw record's layout.
            bool hasBurnIn = value.BurnInEnergyPerDraw is not null || value.BurnInMagPerDraw is not null
                || value.BurnInSumSqClusterSizesPerDraw is not null || value.BurnInMaxClusterSizePerDraw is not null;
            w.Write(hasBurnIn);
            if (hasBurnIn)
                WriteRecord(w, "Burn-in trace",
                    value.BurnInEnergyPerDraw, value.BurnInMagPerDraw,
                    value.BurnInSumSqClusterSizesPerDraw, value.BurnInMaxClusterSizePerDraw);

            // Exact cluster-size histogram (v11): presence flag + its N + 1 counts.
            bool hasExact = value.ExactClusterSizeHistogram is not null;
            w.Write(hasExact);
            if (hasExact)
            {
                w.Write(value.ExactClusterSizeHistogram!.Length);
                foreach (long x in value.ExactClusterSizeHistogram) w.Write(x);
            }

            w.Write(value.TimestampUtc.Ticks);
        }

        /// <summary>
        /// One per-draw record: its length, then E, m, Σ|c|² and c_max, one series after another. The four series
        /// must be present with one length.
        /// </summary>
        private static void WriteRecord(
            BinaryWriter w, string what, double[]? energy, double[]? mag, long[]? sumSq, int[]? maxSize)
        {
            if (energy is null || mag is null || sumSq is null || maxSize is null
                || mag.Length != energy.Length || sumSq.Length != energy.Length || maxSize.Length != energy.Length)
                throw new ArgumentException(
                    $"{what} series are not all present with one length (E {energy?.Length.ToString() ?? "null"}, " +
                    $"m {mag?.Length.ToString() ?? "null"}, Σ|c|² {sumSq?.Length.ToString() ?? "null"}, " +
                    $"c_max {maxSize?.Length.ToString() ?? "null"}).");

            w.Write(energy.Length);
            foreach (double x in energy) w.Write(x);
            foreach (double x in mag) w.Write(x);
            foreach (long x in sumSq) w.Write(x);
            foreach (int x in maxSize) w.Write(x);
        }

        private static (double[] Energy, double[] Mag, long[] SumSq, int[] MaxSize) ReadRecord(BinaryReader r)
        {
            int records = r.ReadInt32();
            var energy = new double[records];
            for (int i = 0; i < records; i++) energy[i] = r.ReadDouble();
            var mag = new double[records];
            for (int i = 0; i < records; i++) mag[i] = r.ReadDouble();
            var sumSq = new long[records];
            for (int i = 0; i < records; i++) sumSq[i] = r.ReadInt64();
            var maxSize = new int[records];
            for (int i = 0; i < records; i++) maxSize[i] = r.ReadInt32();
            return (energy, mag, sumSq, maxSize);
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
            var hist = new long[histBins];
            for (int i = 0; i < histBins; i++) hist[i] = r.ReadInt64();

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

            // Per-draw record (v9).
            var record = ReadRecord(r);

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

            // Burn-in trace (v10).
            (double[] Energy, double[] Mag, long[] SumSq, int[] MaxSize)? burnIn = null;
            if (r.ReadBoolean())
                burnIn = ReadRecord(r);

            // Exact cluster-size histogram (v11).
            long[]? exact = null;
            bool hasExact = r.ReadBoolean();
            if (hasExact)
            {
                int len = r.ReadInt32();
                exact = new long[len];
                for (int i = 0; i < len; i++) exact[i] = r.ReadInt64();
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
                ExactClusterSizeHistogram = exact,
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
                EnergyPerDraw            = record.Energy,
                MagPerDraw               = record.Mag,
                SumSqClusterSizesPerDraw = record.SumSq,
                MaxClusterSizePerDraw    = record.MaxSize,
                BurnInEnergyPerDraw            = burnIn?.Energy,
                BurnInMagPerDraw               = burnIn?.Mag,
                BurnInSumSqClusterSizesPerDraw = burnIn?.SumSq,
                BurnInMaxClusterSizePerDraw    = burnIn?.MaxSize,
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
