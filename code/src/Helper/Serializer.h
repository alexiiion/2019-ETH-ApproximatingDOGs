#pragma once

#include <string>
class DataModel;
struct ViewModel;

void serialize(const std::string& filename, const DataModel& data_model);
void deserialize(const std::string& filename, DataModel& data_model, ViewModel& view_model, bool is_absolute_path = false);
