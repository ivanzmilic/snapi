#ifndef __CMDCFG_H__ // __CMDCFG_H__
#define __CMDCFG_H__

#define CMDOPTS {\
                 "-h|--help            os\"<none>\"                      h\"Usage: momfbd_slave [options]. Options are:\n\n    -p --port\n    -m --master  <machine>\n    -h --help [option]\n    -q --quiet\n    -v --verbose -V\n\nUse -h <option> for more details.\n\"",\
                 "-m|--master          s\"localhost\"                    h\"-m --master: name of machine with master on it [localhost]\"",\
                 "-p|--port            i\"5100\"                         h\"-p --port: port on which master is listening [5100]\"",\
                 "-n|--nthread         i\"1\"                            h\"-n --nthread: number of threads to start [1]\"",\
                 "-to|--timeout        d\"10.0\"                         h\"-to --timeout: time before retrying to connect to master [10.0]\"",\
                 "-t|--timestamp                                         h\"-t --timestamp: print time for all messages\"",\
                 "-q|--quiet                                             h\"-q --quiet: decrease verbosity\"",\
                 "-v|--verbose                                           h\"-v --verbose: increase verbosity\"",\
                 "-V|--version                                           h\"-V: display the version of momfbd\"",\
                 0}

#define CMDMAPS {0}

#endif                // __CMDCFG_H__
