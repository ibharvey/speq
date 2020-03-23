
// CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(
//     std::filesystem::path,
//     cereal::specialization::non_member_load_save_minimal
// );

// template<class Archive>
// void serialize(Archive& archive, std::filesystem::path a)
// {
//     archive(a.string());
// }