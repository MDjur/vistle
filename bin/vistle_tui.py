#! /bin/sh

import argparse
import vistle

console_functions = [
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
    vistle.Status
]

parser = argparse.ArgumentParser(description="Vistle console")
parser.add_argument("host", nargs=1, help="connection host.")
parser.add_argument("port", nargs=1, help="connection port.")
for f in console_functions:
    name = f.__name__
    doc = f.__doc__
    # tmp = [name[0]]
    # shortcut = "".join(tmp + [c for c in name if c.isupper()])
    # parser.add_argument("-" + shortcut,"--" + name, nargs="+", help=doc )
    parser.add_argument("--" + name, nargs="+", metavar="name", help=doc )

args = parser.parse_args()

PORT = args.port[0]
HOST = args.host[0]

vistle._vistle.sessionConnect(None, HOST, PORT)

# for arg in args:
#     print(arg.__name__)
