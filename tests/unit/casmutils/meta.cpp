#include "../../autotools.hh"
#include <filesystem>
#include <gtest/gtest.h>

TEST(AutoToolsTest, InputFilePath)
{
    casmutils::fs::path input_files_path(casmutils::autotools::input_filesdir);
    std::cout << input_files_path << std::endl;
    EXPECT_TRUE(casmutils::fs::exists(input_files_path / "meta.txt"));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
