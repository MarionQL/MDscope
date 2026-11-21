"""Expected, user-correctable failures reported by the CLI without a traceback.

Raise these only at input boundaries (config loading, topology/selection
preflight, CLI arguments). Each subclass keeps its builtin base so library
callers that catch ``ValueError``/``FileNotFoundError``/``IndexError`` still work.
"""


class UserInputError(Exception):
    """Base marker; ``mdscope`` prints ``mdscope: error: ...`` and exits 2."""


class ConfigError(UserInputError, ValueError):
    """Invalid configuration file, option, or topology selection."""


class InputPathError(UserInputError, FileNotFoundError):
    """A configuration or configured input file does not exist."""


class ReplicateIndexError(UserInputError, IndexError):
    """A replicate index is outside the configured system x replicate list."""
