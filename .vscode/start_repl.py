# This script starts an interactive Python REPL with the project in the Python path
import code
import sys
import os
import atexit

# Set up history file - create a custom history file in the project directory
# to avoid issues with the system-wide history file
history_dir = os.path.join(os.path.dirname(__file__), '.history')
os.makedirs(history_dir, exist_ok=True)
histfile = os.path.join(history_dir, 'python_history')

# Try to use readline for history and tab completion
try:
    import readline
    import rlcompleter

    # Set up history
    try:
        if os.path.exists(histfile):
            readline.read_history_file(histfile)
        # Default history len is -1 (infinite), which may grow unruly
        readline.set_history_length(1000)
        # Save history on exit
        atexit.register(readline.write_history_file, histfile)
    except (FileNotFoundError, OSError, IOError) as e:
        print(f"Warning: Could not set up history file: {e}")

    # Enable tab completion
    readline.parse_and_bind("tab: complete")

    print("Command history and tab completion enabled.")
except ImportError:
    print("Module readline not available - history and tab completion disabled.")

# Add the project root to the Python path if it's not already there
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# Print some helpful information
print("Interactive Python REPL with debugging enabled")
print(f"Project root: {project_root}")
print("You can import project modules directly, e.g.:")
print("from syotools.models.telescope import Telescope")
print("\nBreakpoints in your code will be hit when you call functions.")

# Start the interactive console with a more useful set of local variables
namespace = globals().copy()
namespace.update(locals())
code.interact(local=namespace)