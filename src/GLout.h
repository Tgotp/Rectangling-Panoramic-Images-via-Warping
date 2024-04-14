
#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include <opencv2/opencv.hpp>
#include "shader_s.h"
#include "stb_image.h"
#include <iostream>
using namespace std;
using namespace cv;

GLuint matToTexture(const cv::Mat &mat, GLenum minFilter, GLenum magFilter, GLenum wrapFilter) ;
void GLout(Mat img,vector<vector<Point> > mesh);
int GLinit();