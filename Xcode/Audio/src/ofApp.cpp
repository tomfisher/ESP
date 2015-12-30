#include "ofApp.h"

#include <algorithm>

#include "user.h"

//--------------------------------------------------------------
void ofApp::setup() {
    buffer_size_ = 256;
    is_recording_ = false;

    istream_.reset(new AudioStream());
    istream_->onDataReadyEvent(this, &ofApp::onDataIn);

    // setupPipeline is a user-defined function.
    pipeline_ = setupPipeline();

    plot_inputs_.setup(buffer_size_, 1, "Input");
    plot_inputs_.setDrawGrid(true);
    plot_inputs_.setDrawInfoText(true);

    size_t num_pre_processing = pipeline_.getNumPreProcessingModules();
    if (num_pre_processing > 0) {
        PreProcessing* pp = pipeline_.getPreProcessingModule(0);
        plot_pre_processed_.setup(buffer_size_, pp->getNumOutputDimensions(),
                                  "PreProcessing");
        plot_pre_processed_.setDrawGrid(true);
        plot_pre_processed_.setDrawInfoText(true);
    }

    size_t num_feature_modules = pipeline_.getNumFeatureExtractionModules();
    if (num_feature_modules > 0) {
        FeatureExtraction* fe = pipeline_.getFeatureExtractionModule(1);
        for (int i = 0; i < fe->getNumOutputDimensions(); i++) {
            GRT::ofxGrtTimeseriesPlot plot;
            plot.setup(buffer_size_, 1, "Feature");
            plot.setDrawGrid(true);
            plot.setDrawInfoText(true);
            plot_features_.push_back(plot);
        }
    }

    training_data_.setNumDimensions(1);
    training_data_.setDatasetName("Audio");
    training_data_.setInfoText("This data contains audio data");

    ofBackground(54, 54, 54);
}

//--------------------------------------------------------------
void ofApp::update() {
    for (int i = 0; i < input_data_.size(); i++){
        vector<double> data_point = {
            input_data_[i]
        };

        plot_inputs_.update(data_point);

        if (istream_->hasStarted()) {
            if (!pipeline_.preProcessData(data_point)) {
                ofLog(OF_LOG_ERROR) << "ERROR: Failed to compute features!";
            }

            // vector<double> pre_processed_data = pipeline_.getPreProcessedData();
            // plot_pre_processed_.update(pre_processed_data);

            // // The feature vector might be of arbitrary size depending
            // // on the feature selected. But each one could simply be a
            // // time-series.
            // vector<double> feature = pipeline_.getFeatureExtractionData();

            // for (int i = 0; i < feature.size(); i++) {
            //     vector<double> v = { feature[i] };
            //     plot_features_[i].update(v);
            // }
        }

        if (is_recording_) {
            sample_data_.push_back(data_point);
        }
    }
}

//--------------------------------------------------------------
void ofApp::draw() {
    ofSetColor(255);

    int plotX = 10;
    int plotY = 30;
    int plotW = ofGetWidth() - plotX * 2;
    int plotH = 200;
    int margin = 10;

    ofPushStyle();
    ofPushMatrix();
    {
        ofDrawBitmapString("Input:", plotX, plotY - margin);
        plot_inputs_.draw(plotX, plotY, plotW, plotH);
        plotY += plotH + 3 * margin;
    }
    ofPopStyle();
    ofPopMatrix();

    ofPushStyle();
    ofPushMatrix();
    {
        ofDrawBitmapString("PreProcessed:", plotX, plotY - margin);
        plot_pre_processed_.draw(plotX, plotY, plotW, plotH);
        plotY += plotH + 3 * margin;
    }
    ofPopStyle();
    ofPopMatrix();

    ofPushStyle();
    ofPushMatrix();
    ofDrawBitmapString("Feature:", plotX, plotY - margin);

    // int width = plotW / plot_features_.size();
    // for (int i = 0; i < plot_features_.size(); i++) {
    //     plot_features_[i].draw(plotX + i * width, plotY, width, plotH);
    // }
    // plotY += plotH + 3 * margin;

    ofPopStyle();
    ofPopMatrix();
}

void ofApp::exit() {
    if (training_thread_.joinable()) {
        training_thread_.join();
    }
    istream_->stop();
}

void ofApp::onDataIn(vector<double> input) {
    std::lock_guard<std::mutex> guard(input_data_mutex_);
    input_data_ = input;
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
    if (key >= '0' && key <= '9') {
        is_recording_ = true;
        label_ = key - '0';
        sample_data_.clear();
    }

    switch (key) {
        case 't': {
            // If prior training has not finished, we wait.
            if (training_thread_.joinable()) {
                training_thread_.join();
            }

            GRT::TimeSeriesClassificationData data_copy = training_data_;
            auto training_func = [this, data_copy]() mutable {
                ofLog() << "Training started";
                if (pipeline_.train(data_copy)) {
                    ofLog() << "Training is successful";
                } else {
                    ofLog() << "Failed to train the model";
                }
            };
            training_thread_ = std::thread(training_func);
            // training_func();
            break;
        }

        case 's':
            istream_->start();
            break;
        case 'e':
            istream_->stop();
            input_data_.clear();
            break;
        case 'p':
            sample_data_.clear();
            is_recording_ = true;
            break;
    }
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key) {
    is_recording_ = false;
    if (key >= '0' && key <= '9') {
        training_data_.addSample(label_, sample_data_);
    } else if (key == 'p') {
        pipeline_.predict(sample_data_);
        ofLog() << pipeline_.getPredictedClassLabel();
    }
}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ) {

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button) {

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button) {

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button) {

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h) {
}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg) {

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo) {

}