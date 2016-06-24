#include <assert.h>

#include <fstream>
#include <sstream>
#include <string>

#include "GRT/GRT.h"
#include "MFCC.h"

const uint32_t kAudioSampleNum = 42741;
const uint32_t kMFCCFeatureDim = 13;
const uint32_t kMFCCFrameCount = 265;
const uint32_t kSampleRate = 16000;

vector<string> split(string str, char delimiter);
GRT::MatrixDouble readCSVToMatrix(const string& filename, uint32_t row, uint32_t col);

int main(int argc, const char* argv[]) {
    if (argc != 2) {
        std::cerr << "USAGE: " << std::endl;
        std::cerr << "  " << argv[0] << " <data directory>" << std::endl;
        exit(1);
    }

    std::string dir = argv[1];

    GRT::MatrixDouble fft = readCSVToMatrix(dir + "fft.csv", 257, 264);
    GRT::MatrixDouble fbe = readCSVToMatrix(dir + "fbe.csv", 20, 264);

    // 3. Construct the simple pipeline
    return 0;
}

// FFT MAG: 257, 264
// FBE: 20, 264
// CC: 13, 264

#define EXPECT_EQ(a, b) \
    if ((a) != (b)) {                           \
        std::cerr << (a) << " != " << (b) << std::endl; \
        assert(false);                                  \
    }

GRT::MatrixDouble readCSVToMatrix(const string& filename, uint32_t row, uint32_t col) {
    std::cerr << "Reading file: " << filename << std::endl;
    GRT::MatrixDouble mat(row, col);
    std::ifstream file(filename);
    uint32_t r = 0, c = 0;
    string line;
    while (std::getline(file, line)) {
        c = 0;
        // parse csv file
        vector<string> elems = split(line, ',');
        for (const string& elem : elems) {
            mat[r][c++] = std::stod(elem);
        }
        EXPECT_EQ(c, col);
        r++;
    }

    EXPECT_EQ(r, row);
    return mat;
}

vector<string> split(string str, char delimiter) {
    vector<string> internal;
    stringstream ss(str);
    string item;

    while(getline(ss, item, delimiter)) {
        internal.push_back(item);
    }

    return internal;
}
