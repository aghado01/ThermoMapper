using System.Linq;
using Archivory.Tabular;
using Xunit;

namespace VizCore.Tests;

public sealed class TabularProjectionTests
{
    [Fact]
    public void ProjectToTabular_CreatesTableFromGenericSequence()
    {
        var items = new[]
        {
            new { Id = 1, Name = "alpha" },
            new { Id = 2, Name = "beta" },
        };

        var projection = items.ProjectToTabular(
            tableName: "items",
            columns: new[] { "Id", "Name" },
            item => item.Id,
            item => item.Name);

        string csv = projection.ToCsv();

        Assert.Equal("items", projection.TableName);
        Assert.Equal(2, projection.Columns.Count);
        Assert.Contains("Id,Name", csv);
        Assert.Contains("1,alpha", csv);
        Assert.Contains("2,beta", csv);
    }
}
