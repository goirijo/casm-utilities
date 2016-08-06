#ifndef CASM_UTILS_HANDLERS
#define CASM_UTILS_HANDLERS

#include <casm/completer/Handlers.hh>

/**
 * An extension of CASM::OptionHandlerBase, repurposed
 * for bash completion functionality of casm-utilities
 * instead of casm.
 */

namespace casmUtilitiesCompletion
{
    class SplashOption : public CASM::Completer::OptionHandlerBase
    {
        public:

            SplashOption();

            int number() const;

        private:

            void initialize() override;

            int m_number;
    };

    class StructureOption : public CASM::Completer::OptionHandlerBase
    {
        public:

            using CASM::Completer::OptionHandlerBase::output_path;

            StructureOption();

        private:

            void initialize() override;
    };

}

#endif

