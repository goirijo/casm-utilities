#ifndef HANDLERS
#define HANDLERS

#include <casm/completer/Handlers.hh>

/**
 * An extension of CASM::OptionHandlerBase, repurposed
 * for bash completion functionality of casm-utilities
 * instead of casm.
 */

namespace casmUtilities
{
    class structureOption : public CASM::Completer::OptionHandlerBase
    {
        public:

            using CASM::Completer::OptionHandlerBase::output_path;

        private:

            void initialize() override;
    };
}

#endif

