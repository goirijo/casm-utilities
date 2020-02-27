#ifndef CASM_UTILS_EXCPEPTIONS_HH
#define CASM_UTILS_EXCPEPTIONS_HH

#include <casmutils/definitions.hpp>
#include <stdexcept>

namespace except
{

class NotImplemented : public std::runtime_error
{
public:
    NotImplemented() : std::runtime_error("Routine has not been implemented") {}
};

class IncompatibleCoordinate : public std::runtime_error
{
public:
    IncompatibleCoordinate() : std::runtime_error("The provided coordinate is incompatible or could not be mapped.") {}
};

class BasisMismatch : public std::runtime_error
{
public:
    BasisMismatch() : std::runtime_error("The basis between two structures is not compatible") {}

private:
};

class BadCoordMode : public std::runtime_error
{
public:
    BadCoordMode() : std::runtime_error("Invalid mode specified, options are 'CART' and 'FRAC'.") {}

private:
};
class BadPath : public std::runtime_error
{
public:
    BadPath(const casmutils::fs::path& path) : std::runtime_error("Invalid path specified, Received: " + path.string())
    {
    }

private:
};
class OverwriteException : public std::runtime_error
{
public:
    OverwriteException(std::string target_path)
        : std::runtime_error("The path to " + target_path +
                             " already exists! Exception thrown to avoid overwritting files")
    {
    }
};

class UserInputMangle : public std::runtime_error
{
public:
    UserInputMangle(const std::string& init_message) : std::runtime_error(init_message) {}
};
} // namespace except
#endif
