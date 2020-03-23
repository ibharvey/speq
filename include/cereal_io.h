#pragma once

#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/chrono.hpp>
#include <cereal/types/string.hpp>


namespace cereal
{
    template<class Archive>
    void CEREAL_LOAD_MINIMAL_FUNCTION_NAME(
        const Archive&, std::filesystem::path& out, const std::string& in
    )
    {
        out = in;
    }

    template<class Archive>
    std::string CEREAL_SAVE_MINIMAL_FUNCTION_NAME(
        const Archive& ar, const std::filesystem::path& p
    )
    {
        return p.string();
    }
}