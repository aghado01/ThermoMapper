using System.Collections.Generic;

namespace UserRepl;

public sealed record SpcGeneratorInfo(
    string GeneratorName,
    string TypeName,
    string? Description,
    IReadOnlyList<SpcGeneratorParameter> Parameters);

public sealed record SpcGeneratorParameter(
    string Name,
    string TypeName,
    bool IsOptional,
    object? DefaultValue,
    string? Description);
