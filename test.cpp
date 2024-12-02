#include <iostream>
#include <string>
#include <regex>

int main() {
    std::string str = "Tet cells : 0";
    std::regex pattern(R"((\d+))"); // 匹配一个或多个数字
    std::smatch matches;

    if (std::regex_search(str, matches, pattern)) {
        int number = std::stoi(matches[1].str());
        std::cout << "number is: " << number << std::endl;
    } else {
        std::cout << "can't find!" << std::endl;
    }

    return 0;
}