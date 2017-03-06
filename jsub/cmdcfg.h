#ifndef __CMDCFG_H__ // __CMDCFG_H__
#define __CMDCFG_H__

#define CMDOPTS {\
                 "-cfg|--config        s\"mfbd.cfg\"                     h\"-cfg --config <filename>: specify an alternative configuration file from the default \"defaults.dat\".\nThis file contains command line options as they should be specified on the command line. It is read before processing the command line and can be used to change the built-in default values of various program parameters.\"",\
                 "-f                                                     h\"-f: overwrite output file if exists\"",\
                 "-h|--help            os\"<none>\"                      h\"Usage: jsub [options]. Options are:\n\n    -cfg --config <filename>\n    -f\n    -h --help [option]\n    -cl --level <number>\n    -lg <name>\n    -m <master>\n    -n --number <number>\n    -name <job name>\n    -p <port>\n    -pri <priority>\n    -q --quiet\n    -v --verbose -V\n\nUse -h <option> for more details.\n\"",\
                 "-cl|--level          i\"3\"                            h\"-cl --level <number>: configuration data compression level (0..9, higher=better, lower=faster) [3]\"",\
                 "-cc|--chunk_compress i\"1\"                            h\"-cc --chunk_compress <number>: chunk data compression level (0..9, higher=better, lower=faster) [1]\"",\
                 "-lg                  s\"<none>\"                       h\"-lg <name>: log file name\"",\
                 "-m|--master          s\"localhost\"                    h\"-m --master: name of machine with master on it [localhost]\"",\
                 "-mt|--master_threads i\"1\"                            h\"-mt --master_threads: number of threads to use when preprocessing a job on the master [1]\"",\
                 "-ns|--num_slaves     i\"1\"                            h\"-ns --num_slaves: number of slaves to use for this job [1]\"",\
                 "-n|--numbers         i\"<none>\"                       h\"-n --number <number>: Include image sequence number <number>. There is no default.\"",\
                 "-name                s\"<none>\"                       h\"-name <name>: give job a name.\"",\
                 "-o|--output-file     s\"<none>\"                       h\"-p --port: port on which master is listening [5100]\"",\
                 "-p|--port            i\"5100\"                         h\"-o --output-file:  comma separated list of output file base names. File names are applied to the objects in the order they are found in the config file. If insufficient names are provided, a default name will be created for the remaining objects. Excess names will generate a warning but are ignored otherwise. The names are base names only, an appropriate suffix will be attached (.fits/.ana).\"",\
                 "-pri                 i\"0\"                            h\"-pri: priority, a number between -1000 and 1000 [0]\"",\
                 "-q|--quiet                                             h\"-q --quiet: decrease verbosity\"",\
                 "-s|--swap                                              h\"-s --swap: swap mode: write compressed data to swap file instead of keeping it in memory (useful for large problems).\"",\
                 "-st|--slave_threads  i\"0\"                            h\"-st --slave_threads: maximum number of threads to use when processing a chunk on the slaves [maximum allowed by slave]\"",\
                 "-S                   os\"1.0,0.0,0.0,0.0\"             h\"-S: Stokes flatfielding values.\"",\
                 "-t|--timestamp                                         h\"-t --timestamp: print time for all messages\"",\
                 "-v|--verbose                                           h\"-v --verbose: increase verbosity\"",\
                 "-V|--version                                           h\"-V: display the version of momfbd\"",\
                 0}

#define CMDMAPS {0}

#endif                // __CMDCFG_H__
