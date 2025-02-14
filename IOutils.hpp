#include <fstream>
#include <cstring>
#include <vector>
#include <iostream>
#include <cmath>

void load_fasta(std::string filename, std::vector<char> &Text, uint64_t &sum, uint64_t &ns)
{
    // open the input stream
    std::ifstream input(filename);
    if (!input.good())
    {
        std::cerr << "Error opening " << filename << ". exiting..." << std::endl;
        exit(-1);
    }
    input.seekg(0, std::ios::end);
    Text.resize(input.tellg());
    input.seekg(0, std::ios::beg);

    std::string line, DNA_sequence;
    sum = 0, ns = 0;
    while (std::getline(input, line))
    {
        // skip empty lines
        if (line.empty())
        {
            continue;
        }
        // header of a new sequence
        if (line[0] == '>')
        {
            ns++; // increase sequence count
            size_t seqlen = DNA_sequence.size();
            // insert previous DNA sequence
            if (seqlen > 0)
            {
                // 1. 反转 DNA_sequence
                std::reverse(DNA_sequence.begin(), DNA_sequence.end());

                // 2. 将翻转后的序列写入 Text
                memcpy(&Text[sum], &DNA_sequence[0], seqlen);
                // increase current text size
                sum += seqlen;
                // add separator
                Text[sum++] = 1;
            }
            // the current DNA sequence
            DNA_sequence.clear();
        }
        else
        {
            DNA_sequence += line;
        }
    }
    memcpy(&Text[sum], &DNA_sequence[0], DNA_sequence.size());
    sum += DNA_sequence.size();
    Text[sum++] = 1;
    Text.resize(sum);
    Text.shrink_to_fit();
    // close stream
    input.close();
}

void writeTextToFile(const std::string &filename, const std::vector<char> &Text)
{
    // 创建一个输出文件流，打开文件以写入二进制数据
    std::ofstream outfile(filename, std::ios::binary);

    // 检查文件是否成功打开
    if (!outfile)
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // 写入数据
    if (!Text.empty())
    {
        outfile.write(reinterpret_cast<const char *>(Text.data()), Text.size());

        // 检查写入是否成功
        if (!outfile)
        {
            std::cerr << "Error: Failed to write data to " << filename << "." << std::endl;
        }
    }
    else
    {
        std::cerr << "Warning: Text vector is empty. Nothing was written to " << filename << "." << std::endl;
    }

    // 关闭文件流
    outfile.close();
}