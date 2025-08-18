import os, glob, re, threading, queue
from collections import deque
import subprocess as sub
import numpy as np
from util.qol_utils.progress_bar import printProgressBar, printProgressWithOutput
from typing import Optional

class DelphesSetup:
    """
    A class for running downloading/building (if necessary) the
    Delphes fast detector simulation framework.
    """
    def __init__(self,delphes_dir:Optional[str]=None):
        self.SetDirectory(delphes_dir)
        self.executable = None
        self.prefix = self.__class__.__name__
        self.expected_stdout = 281
        self.build_dir = None

    def SetDirectory(self,delphes_dir:Optional[str]=None):
        self.delphes_dir = delphes_dir
        if(self.delphes_dir is None):
            self.delphes_dir = os.path.dirname(os.path.abspath(__file__)) + '/../../external/delphes'
        self.delphes_dir = os.path.normpath(self.delphes_dir)

        # Make the Delphes dir if it does not exist.
        os.makedirs(self.delphes_dir,exist_ok=True)

        # Files that are used for logging progress with downloading/building Delphes.
        self.logfile = '{}/log.stdout'.format(self.delphes_dir)
        self.errfile = '{}/log.stderr'.format(self.delphes_dir)

    def GetDirectory(self)->str:
        return self.delphes_dir

    def GetExecutable(self)->str:
        return self.executable

    def Prepare(self, j:int=4, force:bool=False, verbose:bool=False):
        """
        Checks if Delphes is built. If it's not,
        it will build Delphes.
        (Historically, this class handled downloading code as well,
         but we have shifted to using git submodules for better versioning.)
        """
        # Check if Delphes is already built at destination.
        if(not force):
            ex1 = glob.glob('{}/**/DelphesHepMC3'.format(self.delphes_dir),recursive=True)
            if(len(ex1) > 0):
                self.executable = ex1[0]
                if(verbose): self._print('Found DelphesHepMC3 @ {}'.format(self.executable['hepmc']))
                return

        self.Build(j=j)
        self.executable = glob.glob('{}/**/DelphesHepMC3'.format(self.delphes_dir),recursive=True)[0]
        return

    def Build(self, j:int=4):
        """
        Builds Delphes. The official documentation suggests simply running "make", which will use the existing Makefile
        that ships with Delphes. However, we're going to invoke CMake since Delphes also ships a CMakeLists.txt file,
        and it looks like this will automatically build the display library (that we may also want to leverage).
        """

        if(self.build_dir is None):
            self.build_dir = self.delphes_dir

        self._print('Configuring Delphes.')
        command = ['cmake','./']
        self.run_command_with_progress(command,cwd=self.build_dir, prefix='Configuring Delphes:',output_length=18)

        self._print('Building Delphes @ {} .'.format(self.build_dir))
        command = ['cmake', '--build','.', '-j{}'.format(j)]
        self.run_command_with_progress(command,cwd=self.build_dir, prefix='Building Delphes:')#,output_length=self.expected_stdout)
        return

    def _print(self,val:str):
        print('{}: {}'.format(self.prefix,val))

    def run_command_with_progress(self, command, prefix, cwd=None,
                                show_output_lines=10, output_width=80, output_length=None):
            """
            Run command with progress bar and scrolling output display

            Args:
                command: Command to run
                cwd: Current working directory for command
                prefix: Progress bar prefix
                show_output_lines: Number of recent output lines to show (0 to disable)
                output_width: Width of output display area
                output_length: Optional expected length of output; to use if output has no progress tracking (e.g. "[1%] ...")
            """
            # Pattern to match progress indicators
            progress_pattern = re.compile(r'\[\s*(\d+)%\]')

            # Buffer to store recent output lines
            recent_lines = deque(maxlen=show_output_lines * 2)  # Store more than we show for better filtering

            # Track current progress
            current_progress = 0

            def read_stream(stream, file_handle, line_queue, stream_name):
                """Helper function to read from stdout/stderr in separate threads"""
                try:
                    for line in iter(stream.readline, ''):
                        if line:
                            file_handle.write(line)
                            file_handle.flush()
                            line_queue.put((stream_name, line))
                except Exception as e:
                    line_queue.put(('error', f"Error reading {stream_name}: {str(e)}\n"))
                finally:
                    stream.close()

            with open(self.logfile, 'a') as f, open(self.errfile, 'a') as g:

                # add an empty space -- nice to separate outputs from different commands
                f.write('\n')
                g.write('\n')

                process = sub.Popen(
                    command,
                    cwd=cwd,
                    stdout=sub.PIPE,
                    stderr=sub.PIPE,
                    universal_newlines=True,
                    bufsize=1  # Line buffered
                )

                # Create queues for inter-thread communication
                line_queue = queue.Queue()

                # Start threads to read stdout and stderr
                stdout_thread = threading.Thread(
                    target=read_stream,
                    args=(process.stdout, f, line_queue, 'stdout')
                )
                stderr_thread = threading.Thread(
                    target=read_stream,
                    args=(process.stderr, g, line_queue, 'stderr')
                )

                stdout_thread.daemon = True
                stderr_thread.daemon = True
                stdout_thread.start()
                stderr_thread.start()

                # Initialize display area
                if show_output_lines > 0:
                    # Reserve space for progress bar + output lines
                    for i in range(3 + show_output_lines):
                        print()
                    # Initial display
                    printProgressWithOutput(current_progress, 100, recent_lines,
                                        prefix=prefix, suffix='Complete',
                                        max_lines=show_output_lines, width=output_width)
                else:
                    # Just show progress bar
                    printProgressBar(current_progress, 100, prefix=prefix, suffix='Complete')

                counter = 0
                try:
                    while process.poll() is None or not line_queue.empty():
                        try:
                            # Get line with timeout to check process status
                            stream_name, line = line_queue.get(timeout=0.1)

                            # Add to recent lines buffer
                            recent_lines.append(line)

                            # Update progress -- either checking for progress indicator,
                            # or expecting some fixed length of output.
                            if(output_length is not None and stream_name == 'stdout'):
                                current_progress = (counter / output_length) * 100.
                                current_progress = np.minimum(current_progress, 100.) # for safety

                            # Check for progress indicator in stdout
                            elif stream_name == 'stdout':
                                match = progress_pattern.search(line)
                                if match:
                                    current_progress = int(match.group(1))

                            # Update display (both progress and output)
                            if show_output_lines > 0:
                                printProgressWithOutput(current_progress, 100, recent_lines,
                                                    prefix=prefix, suffix='Complete',
                                                    max_lines=show_output_lines, width=output_width)
                            else:
                                printProgressBar(current_progress, 100, prefix=prefix, suffix='Complete')

                            if(stream_name == 'stdout'): counter += 1
                        except queue.Empty:
                            continue  # Continue checking if process is still running

                except KeyboardInterrupt:
                    process.terminate()
                    raise

                # Wait for threads to complete
                stdout_thread.join(timeout=1.0)
                stderr_thread.join(timeout=1.0)

                # Final progress update
                if show_output_lines > 0:
                    printProgressWithOutput(100, 100, recent_lines,
                                        prefix=prefix, suffix='Complete',
                                        max_lines=show_output_lines, width=output_width)
                    print()  # Move cursor to next line after final display
                else:
                    printProgressBar(100, 100, prefix=prefix, suffix='Complete')
                    print()  # Move to next line

                # Wait for process to complete and check return code
                return_code = process.wait()
                if return_code != 0:
                    raise sub.CalledProcessError(return_code, command)

class DelphesROOTHepMC3Setup(DelphesSetup):
    """
    A class for running building the
    DelphesHepMC3ROOT executable -- our own custom one
    for handling HepMC3/ROOT files.
    """
    def __init__(self):
        self.SetDirectory()
        self.executable = None
        self.prefix = self.__class__.__name__

    def SetDirectory(self,delphes_dir:Optional[str]=None):
        self.delphes_dir = delphes_dir
        if(self.delphes_dir is None):
            self.delphes_dir = os.path.dirname(os.path.abspath(__file__)) + '/../../external/delphes_custom'
        self.delphes_dir = os.path.normpath(self.delphes_dir)

        # # Make the Delphes dir if it does not exist. # Commented out for now, while delphes_custom is part of this package directly
        # os.makedirs(self.delphes_dir,exist_ok=True)

        # Files that are used for logging progress with building DelphesHepMC3ROOT.
        self.logfile = '{}/log.stdout'.format(self.delphes_dir)
        self.errfile = '{}/log.stderr'.format(self.delphes_dir)

    def GetDirectory(self)->str:
        return self.delphes_dir

    def GetExecutable(self)->str:
        return self.executable

    def Prepare(self, j:int=4, force:bool=False, verbose:bool=False):
        """
        Checks if DelphesHepMC3ROOT is built. If it's not,
        it will build it.
        """
        # Check if DelphesHepMC3ROOT is already built at destination.
        if(not force):
            ex1 = glob.glob('{}/**/DelphesHepMC3ROOT'.format(self.delphes_dir),recursive=True)
            if(len(ex1) > 0):
                self.executable = ex1[0]
                if(verbose): self._print('Found DelphesHepMC3ROOT @ {}'.format(self.executable))
                return


        self.Build(j=j)
        self.executable = glob.glob('{}/**/DelphesHepMC3ROOT'.format(self.delphes_dir),recursive=True)[0]
        return

    def Build(self, j:int=4):
        self.build_dir = '{}/{}'.format(self.delphes_dir,'build')
        os.makedirs(self.build_dir,exist_ok=True)
        self._print('Building DelphesHepMC3ROOT @ {} .'.format(self.build_dir))
        command = ['cmake', '../']
        self.run_command_with_progress(command,prefix='Configuring DelphesHepMC3ROOT',cwd=self.build_dir)
        command = ['make', '-j{}'.format(j)]
        self.run_command_with_progress(command,prefix='Building DelphesHepMC3ROOT',cwd=self.build_dir,show_output_lines=5)
        return

    def _print(self,val:str):
        print('{}: {}'.format(self.prefix,val))