import os,glob,pathlib,threading, queue, re
from collections import deque
import subprocess as sub
import ROOT as rt
from util.qol_utils.progress_bar import printProgressBar, printProgressWithOutput

class PileupSetup:
    """
    This class is for setting up the C++/ROOT Pileup code.
    """
    def __init__(self,pileup_dir=None):
        self.SetDirectory(pileup_dir)
        self.executable = '/bin/bash'
        self.lib_extensions = ['so','dylib']
        self.status = False

    def SetDirectory(self,dir):
        if(dir is None):
            dir = os.path.dirname(os.path.abspath(__file__)) + '/../root/pileup'
        self.dir = dir
        self.logfile = '{}/log.stdout'.format(self.dir)
        self.errfile = '{}/log.stderr'.format(self.dir)

    def FullPreparation(self):
        self.status = False
        try: self.Load(quiet=True)
        except:
            self.Build()
            self.Load()
        self.status = True
        return

    def Build(self):
        """
        Builds the Pileup C++/ROOT library.
        This runs a shell script.
        """
        # TODO: Consider some fancy stuff for setting Pythia8 environment variables?
        self.build_script = 'build.sh'

        command = ['./{}'.format(self.build_script)]
        self.run_command_with_progress(command,cwd=self.dir,prefix='Building Pileup:',output_width=80)
        return

    def Load(self,quiet=False):
        # Load our custom ROOT library.
        try:
            a = rt.Pileup
            return
        except:
            pass

        # Note that we *also* need to fetch some include files -- this has something to do with using the ROOT interpreter.
        # Also note that we have multiple library paths -- to allow for Linux/macOS compatibility.
        build_dir = self.dir + '/pileup/build'
        include_dir = self.dir + '/pileup/inc'
        custom_lib_paths = [os.path.realpath('{}/lib/libPileup.{}'.format(build_dir,x)) for x in self.lib_extensions]
        custom_inc_paths = glob.glob(include_dir + '/**/*.h',recursive=True)

        # Check for any of the libraries.
        found_libary = False
        for libpath in custom_lib_paths:
            if(pathlib.Path(libpath).exists()):
                custom_lib_path = libpath
                found_libary = True
                break

        if(not found_libary):
            if(not quiet): print('Error: The Pileup lib has not been built!')
            assert False

        for inc_path in custom_inc_paths:
            command = '#include "{}"'.format(inc_path)
            status = rt.gInterpreter.Declare(command)
            assert status, 'The following header file did not load properly: {}'.format(inc_path)

        status = rt.gSystem.Load(custom_lib_path)
        assert status == 0, 'The Pileup lib did not load properly.'
        return

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
                            if(output_length is not None):
                                current_progress = (counter / output_length) * 100.

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

                            counter += 1
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
