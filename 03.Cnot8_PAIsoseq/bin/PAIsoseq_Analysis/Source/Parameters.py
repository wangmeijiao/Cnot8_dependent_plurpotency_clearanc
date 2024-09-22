import yaml
import os
import sys
from Source import Misc

class Parameters():

    def __init__(self, yamlParameterPath):
        self.parameterPath = yamlParameterPath

        # Check if YAML file exists
        Misc.Misc.checkFile(self.parameterPath, 'Parameters YAML File')

        self.parameters = self.parseYaml(self.parameterPath)

        # regDirs contains all dirs that have been created / checked during pipeline session.
        self.regDirs = dict()

    def parseYaml(self, yamlParameterPath):
        with open(yamlParameterPath, 'r') as yamlStream:
            try:
                parameters = yaml.load(yamlStream)
            except yaml.YAMLError:
                sys.stderr.write("Error Reading Parameter YAML at: %s\n".format(yamlParameterPath))

        return parameters

    def getParameters(self):
        return self.parameters

    def getParametersBysubKey(self, key, subkey):
        if key in self.parameters:
            parameters_dic = self.parameters[key]
            if subkey in parameters_dic:
                return self.parameters[key][subkey]
            else:
                print(subkey + 'not found in Parameter.')
                cmd = "which " + subkey
                path_env = os.popen(cmd).read().strip()
                if path_env == "":
                    print(key + 'not found in $PATH.')
                else:
                    self.parameters[key][subkey] = path_env
                return self.parameters[key][subkey]
        else:
            raise KeyError('No Data found for Key for Parameter.')

    def registerDir(self, key, path):
    # Register dir path for key

        if key in self.regDirs:
            raise KeyError('Key for regDict already registered.')
        else:
            self.regDirs[key] = path

    def getRegDir(self, key):

        if not key in self.regDirs:
            return KeyError('Key not in regDict')
        else:
            return self.regDirs[key]

    def checkKey(self, key):
        if key in self.regDirs:
            return True
        else:
            return False

    def getExpName(self):
        try:
            return self.parameters['experimentName']
        except KeyError:
            raise KeyError("No Experiment Defined In parameters.yaml")
