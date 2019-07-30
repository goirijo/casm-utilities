#ifndef CASM_UTILS_EXCPEPTIONS_HH
#define CASM_UTILS_EXCPEPTIONS_HH

#include <stdexcept>

namespace UtilExcept
{

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
} // namespace UtilExcept
#endif
