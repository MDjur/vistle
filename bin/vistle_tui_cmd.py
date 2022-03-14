#! /bin/sh

from cmd import Cmd
import vistle

class Command:
    def __init__(self, function):
        self.function = function

    def function_wrapper(self, *args, **kwargs):
        return self.function(*args, **kwargs)

    def name(self):
        return self.function.__name__

    def help(self):
        print(self.function.__doc__)

class VistleShell(Cmd):
    prompt = '( vistle ) '
    intro = "Welcome to VistleShell. Type ? to list commands"

    def do_exit(self, inp):
        print("Close")
        return True

    def help_exit(self):
        print('exit the application. Shorthand: x q Ctrl-D.')

    def default(self, inp):
        if inp in ('x', 'q'):
            return self.do_exit(inp)

        print("Default: {}".format(inp))

console_functions = [Command(f) for f in [
    vistle.connect,
    vistle.source,
    vistle.removeHub,
    vistle.spawn,
    vistle.spawnAsync,
    vistle.waitForSpawn,
    vistle.kill,
    vistle.disconnect,
    vistle.compute,
    vistle.interrupt,
    vistle.quit,
    vistle.trace,
    vistle.barrier,
    vistle.requestTunnel,
    vistle.removeTunnel,
    vistle.printInfo,
    vistle.printWarning,
    vistle.printError,
    vistle.setStatus,
    vistle.clearStatus,
    vistle.setLoadedFile,
    vistle.getLoadedFile,
    vistle.setParam,
    vistle.setIntParam,
    vistle.setFloatParam,
    vistle.setStringParam,
    vistle.setVectorParam,
    vistle.setIntVectorParam,
    vistle.applyParameters,
    vistle.getAvailable,
    vistle.getRunning,
    vistle.getBusy,
    vistle.getModuleName,
    vistle.getModuleDescription,
    vistle.waitForHub,
    vistle.waitForHubs,
    vistle.waitForNamedHubs,
    vistle.getMasterHub,
    vistle.getVistleSession,
    vistle.getSessionUrl,
    vistle.getAllHubs,
    vistle.getHub,
    vistle.getPortDescription,
    vistle.getInputPorts,
    vistle.getOutputPorts,
    vistle.getConnections,
    vistle.getParameters,
    vistle.getParameterType,
    vistle.getParameterPresentation,
    vistle.getParameterTooltip,
    vistle.isParameterDefault,
    vistle.getIntParam,
    vistle.getFloatParam,
    vistle.getVectorParam,
    vistle.getIntVectorParam,
    vistle.getStringParam,
    vistle.getEscapedStringParam,
    vistle.Id,
    vistle.Message,
    vistle.StateObserver,
    vistle.Text,
    vistle.Status,
    vistle._vistle.sessionConnect
]]

shell = VistleShell()
for f in console_functions:
    setattr(VistleShell, "do_" + f.name(), f.function_wrapper)
    setattr(VistleShell, "help_" + f.name(), f.help)
shell.cmdloop()
