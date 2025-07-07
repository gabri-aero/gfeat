#ifndef _FILE_READER_HPP_
#define _FILE_READER_HPP_

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

class FileReader {
  public:
    static std::vector<std::vector<std::string>> read(std::string filename,
                                                      char delimiter = ',') {
        std::ifstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Error: file was not opened!" << std::endl;
        }

        std::string line;
        std::vector<std::vector<std::string>> data;
        std::vector<std::string> data_row;

        while (getline(file, line)) {
            std::stringstream ss(line);
            std::string item;
            while (getline(ss, item, delimiter)) {
                if (item.size() > 0)
                    data_row.push_back(item);
            }
            if (data_row.size() > 0)
                data.push_back(data_row);
            data_row.clear();
        }

        return data;
    }
};

#endif // _FILE_READER_HPP_