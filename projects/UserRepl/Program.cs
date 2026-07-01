using UserRepl;

namespace Spc.User;

public static class Program
{
    public static int Main(string[] args) => SubcommandRouter.Run(args);
}
