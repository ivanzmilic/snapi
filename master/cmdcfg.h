#ifndef __CMDCFG_H__ // __CMDCFG_H__
#define __CMDCFG_H__

#define CMDOPTS {\
                 "-g|--grid                                              h\"-g --grid: start GRID slave\"",\
                 "-h|--help            os\"<none>\"                      h\"Usage: manager [options]. Options are:\n\n    -p --port\n    -to  <int>\n    -h --help [option]\n    -q --quiet\n    -v --verbose -V\n\nUse -h <option> for more details.\n\"",\
                 "-p|--port            i\"5100\"                         h\"-p --port: port on which master is listening [5100]\"",\
                 "-to|--timeout        i\"7600\"                         h\"-to --timeout: time (in seconds) to wait before declaring slave dead [7600]\"",\
                 "-t|--timestamp                                         h\"-t --timestamp: print time for all messages\"",\
                 "-q|--quiet                                             h\"-q --quiet: decrease verbosity\"",\
                 "-v|--verbose                                           h\"-v --verbose: increase verbosity\"",\
                 "-V|--version                                           h\"-V: display the version of momfbd\"",\
                 0}

#define CMDMAPS {0}

#endif                // __CMDCFG_H__
