"""Errors shared by biological-resource providers."""


class ProviderError(ValueError):
    """Base class for provider failures shown to the user."""


class ProviderUnavailableError(ProviderError):
    """Raised when a configured resource cannot be opened."""


class ProviderSchemaError(ProviderError):
    """Raised when a resource does not satisfy its documented schema."""


class ResourceNotFoundError(ProviderError):
    """Raised when a requested marker set or reference does not exist."""
