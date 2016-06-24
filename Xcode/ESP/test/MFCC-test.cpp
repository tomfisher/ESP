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

#define EXPECT_CLOSE_VEC(a, b, eps)                         \
    for (auto ai = (a).begin(), bi = (b).begin();           \
         ai != (a).end() && bi != (b).end(); ++ai, ++bi) {  \
        if (abs(*ai - *bi) > eps) {                         \
            std::cerr << "diff(" << *ai << ", " << *bi      \
                      << ") > " << eps << std::endl;        \
            assert(false);                                  \
        }                                                   \
    }

#define EXPECT_EQ(a, b)                                 \
    if ((a) != (b)) {                                   \
        std::cerr << (a) << " != " << (b) << std::endl; \
        assert(false);                                  \
    }

vector<string> split(string str, char delimiter);
GRT::MatrixDouble readCSVToMatrix(const string& filename, uint32_t row, uint32_t col);

int main(int argc, const char* argv[]) {
    if (argc != 2) {
        std::cerr << "USAGE: " << std::endl;
        std::cerr << "  " << argv[0] << " <data directory>" << std::endl;
        exit(1);
    }

    std::string dir = argv[1];

    // FFT MAG: 257, 264
    // FBE: 20, 264
    // CC: 13, 264
    GRT::MatrixDouble fft = readCSVToMatrix(dir + "fft.csv", 257, 264);
    GRT::MatrixDouble fbe = readCSVToMatrix(dir + "fbe.csv", 20, 264);
    GRT::MatrixDouble cc = readCSVToMatrix(dir + "cc.csv", 13, 264);

    GRT::MFCC mfcc(300, 3700, 257, kSampleRate, 20, 13);

    // Step 1. triangle filter
    vector<GRT::MelBank> filters = mfcc.getFilters();

    vector<double> filter = filters[0].getFilter();
    uint32_t M = filters.size();

    uint32_t col = 0;
    GRT::VectorDouble fft_frame = fft.getColVector(col);
    GRT::VectorDouble my_fbe(M);
    for (uint32_t i = 0; i < M; i++) {
        double energy = filters[i].filter(fft_frame);
        if (energy == 0) {
            // Prevent log_energy goes to -inf...
            my_fbe[i] = 0;
        } else {
            my_fbe[i] = log(energy);
        }
    }

    GRT::VectorDouble their_fbe = fbe.getColVector(col);
    EXPECT_CLOSE_VEC(my_fbe, their_fbe, 0.01);

    // Step 2. IDCT
    // TODO

    std::cout << "Passed" << std::endl;
    return 0;
}

GRT::MatrixDouble readCSVToMatrix(const string& filename, uint32_t row, uint32_t col) {
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
