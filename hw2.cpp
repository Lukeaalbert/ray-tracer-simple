/*
  CSCI 420 Computer Graphics, USC
  Assignment 2: Roller Coaster
  C/C++ starter code

  Student username: lpalbert
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <cstring>
#include <memory>
#include "openGLHeader.h"
#include "imageIO.h"
#include "glutHeader.h"
#include "openGLMatrix.h"
#include "pipelineProgram.h"
#include "vbo.h"
#include "vao.h"
#include "ebo.h"

#if defined(WIN32) || defined(_WIN32)
  #ifdef _DEBUG
    #pragma comment(lib, "glew32d.lib")
  #else
    #pragma comment(lib, "glew32.lib")
  #endif
#endif

#if defined(WIN32) || defined(_WIN32)
  char shaderBasePath[1024] = SHADER_BASE_PATH;
#else
  char shaderBasePath[1024] = "../openGLHelper";
#endif

using namespace std;

// *****************************************************************************
// ******************** Global Variables (Starter Code) ************************
// *****************************************************************************

int mousePos[2]; // x,y screen coordinates of the current mouse position

int leftMouseButton = 0; // 1 if pressed, 0 if not 
int middleMouseButton = 0; // 1 if pressed, 0 if not
int rightMouseButton = 0; // 1 if pressed, 0 if not

typedef enum { ROTATE, TRANSLATE, SCALE } CONTROL_STATE;
CONTROL_STATE controlState = ROTATE;

// Transformations of the terrain.
float terrainRotate[3] = { 0.0f, 0.0f, 0.0f }; 
// terrainRotate[0] gives the rotation around x-axis (in degrees)
// terrainRotate[1] gives the rotation around y-axis (in degrees)
// terrainRotate[2] gives the rotation around z-axis (in degrees)
float terrainTranslate[3] = { 0.0f, 0.0f, 0.0f };
float terrainScale[3] = { 1.0f, 1.0f, 1.0f };

// Width and height of the OpenGL window, in pixels.
int windowWidth = 1280;
int windowHeight = 720;
// char windowTitle[512] = "CSCI 420 Homework 2";
char windowTitle[512] = "CSCI 420 Homework 2";

// Number of vertices in the single triangle (starter code).
int numVertices;

// CSCI 420 helper classes.
OpenGLMatrix matrix;
PipelineProgram pipelineProgram;
VBO vboVertices;
VBO vboNormals;
VAO vao;

// Represents one spline control point.
struct Point 
{
  float x, y, z;
};

// Contains the control points of the spline.
struct Spline 
{
  int numControlPoints;
  std::vector<Point> points;
} spline;

// *****************************************************************************
// ******************** Global Variables (New) *********************************
// *****************************************************************************

unsigned int screenshotCount = 0;

// one less than the height of the sky
float lightPosition[] = {0.0f, 29.0f, 0.0f};

// text image file paths
const char groundWallsImageFilePath[] = "imgs/ground.jpg";
const char skyImageFilePath[] = "imgs/sky.jpg";

// texture handels (important)
GLuint groundWallsHandle, skyHandle;

// pipeline program used for groundWalls and sky
PipelineProgram texturePipelineProgram;

// VBOS for groundWalls vertices and text coords
VBO vbogroundWallsVertices, vbogroundWallsTexCoords;

// VBOS for sky vertices and text coords
VBO vboSkyVertices, vboSkyTexCoords;

VAO vaogroundWalls, vaoSky;

// catmull-rom basis matrix (column-major)
const float basisMatrix[16] = {
  -0.5, 1.5, -1.5, 0.5,
  1.0, -2.5, 2.0, -0.5,
  -0.5, 0.0, 0.5, 0.0,
  0.0, 1.0, 0.0, 0.0
};

// points along the spline
// x0, y0, z0, x1, y1, z1, ...
vector<float> points;
vector<float> tangents;

// points along the rail.
// p0, p1, ..., p7 represents one rail chunk. 
vector<float> railPoints;
vector<unsigned int> railIndices;
vector<float> railNormals;
EBO eboRailIndeces;
unsigned int numRailPoints;

// default scalar. used for matmull-rom points and tangents.
const float s = 0.5;
// TODO: change when velocity added! just a constant for now.
int animationSpeed = 8;
// Whether the coaster is moving (animation on)
bool animationOn = false;
// up
float worldUp[3] = {0.0f, 1.0f, 0.0f};
// scales the width of the track rail
const float alpha = 0.1f;

// **************************************************************************
// ****************************** Helper Functions **************************
// **************************************************************************
// Multiply C = A * B, where A is a m x p matrix, and B is a p x n matrix.
// All matrices A, B, C must be pre-allocated (say, using malloc or similar).
// The memory storage for C must *not* overlap in memory with either A or B. 
// That is, you **cannot** do C = A * C, or C = C * B. However, A and B can overlap, and so C = A * A is fine, as long as the memory buffer for A is not overlaping in memory with that of C.
// Very important: All matrices are stored in **column-major** format.
// Example. Suppose 
//      [ 1 8 2 ]
//  A = [ 3 5 7 ]
//      [ 0 2 4 ]
//  Then, the storage in memory is
//   1, 3, 0, 8, 5, 2, 2, 7, 4. 
void MultiplyMatrices(int m, int p, int n, const float * A, const float * B, float * C)
{
  for(int i=0; i<m; i++)
  {
    for(int j=0; j<n; j++)
    {
      float entry = 0.0;
      for(int k=0; k<p; k++)
        entry += A[k * m + i] * B[j * p + k];
      C[m * j + i] = entry;
    }
  }
}

float vecLen(const float v[3]) {
  return sqrtf(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}
void normalize(float v[3]) {
  float L = vecLen(v);
  if (L > 1e-8f) {
    v[0] /= L;
    v[1] /= L;
    v[2] /= L;
    }
}
void cross(const float a[3], const float b[3], float out[3]) {
  out[0] = a[1]*b[2] - a[2]*b[1];
  out[1] = a[2]*b[0] - a[0]*b[2];
  out[2] = a[0]*b[1] - a[1]*b[0];
}

void drawTriangle(unsigned int A, unsigned int B, 
  unsigned int C, vector<unsigned int>& indeces) {
    indeces.push_back(A);
    indeces.push_back(B);
    indeces.push_back(C);
}

void drawSquare(unsigned int bottomRightIdx,  unsigned int bottomLeftIdx,
  unsigned int topRightIdx,  unsigned int topLeftIdx,
  vector<unsigned int>& indeces)
{
  drawTriangle(bottomRightIdx, topRightIdx, bottomLeftIdx, indeces);
  drawTriangle(topLeftIdx, bottomLeftIdx, topRightIdx, indeces);
}

// **************************************************************************
// ****************************** Coaster Struct ****************************
// **************************************************************************

struct Coaster {
    // the index into "points". where are on the spline.
    // always a multiple of 3.
    unsigned int coasterPoint = 3;
    float u = 0.0f;
    // move numSteps steps forward in the animation.
    void stepForward(unsigned int numSteps = 1) {
      for (int i = 0; i < numSteps; ++i) {
        if (u >= 1.0f) u = 0.0f;
        else u+=0.001f;
        coasterPoint += 3;
      }
    }
    // get a point on the coaster.
    // offset = the amount we want to look forwards or backwards
    // in "points".
    void getCoasterPoint(float out[3], int offset = 0) {
      out[0] = points[coasterPoint + 3*offset];
      out[1] = points[coasterPoint + 1 + 3*offset];
      out[2] = points[coasterPoint + 2 + 3*offset];
    }
    void lookAtOnCoaster() {
      // set camera to be on the current point on the track, pointing
      // towards the tangent.
      float prevPoint[3], point[3], nextPoint[3], nextNextPoint[3];
      getCoasterPoint(prevPoint, -1);
      getCoasterPoint(point);
      getCoasterPoint(nextPoint, 1);
      getCoasterPoint(nextNextPoint, 2);
      float controlMatrix[12] = {
        prevPoint[0], point[0], nextPoint[0],  nextNextPoint[0],
        prevPoint[1],  point[1], nextPoint[1],  nextNextPoint[1],
        prevPoint[2], point[2], nextPoint[2],  nextNextPoint[2],
      };
      // get tangent
      float tangent[3] = {
        tangents[coasterPoint],
        tangents[coasterPoint + 1],
        tangents[coasterPoint + 2]
      };

      // create lookat vec
      float lookAt[3] = {
      point[0] + 2.0f * tangent[0],
      point[1] + 2.0f * tangent[1],
      point[2] + 2.0f * tangent[2]
      };

      // stable up vector
      float right[3], up[3];
      cross(tangent, worldUp, right);
      normalize(right);
      cross(right, tangent, up);
      normalize(up);

      // set lookat
      matrix.LookAt(
        point[0], point[1]+1.0f, point[2],
        lookAt[0], lookAt[1], lookAt[2],
        up[0], up[1], up[2]
      );

      // send camera position to shaders so we can do specular lighting 
      pipelineProgram.Bind();
      pipelineProgram.SetUniformVariableVec3f("cameraPosition", point[0], point[1]+1.0f, point[2]);
    }
} coaster;

// **************************************************************************
// ****************************** Main Program Logic ************************
// **************************************************************************

// Write a screenshot to the specified filename.
void saveScreenshot(unsigned int screenshotNumber)
{
  std::unique_ptr<unsigned char[]> screenshotData = std::make_unique<unsigned char[]>(windowWidth * windowHeight * 3);
  glReadPixels(0, 0, windowWidth, windowHeight, GL_RGB, GL_UNSIGNED_BYTE, screenshotData.get());

  ImageIO screenshotImg(windowWidth, windowHeight, 3, screenshotData.get());

  // get file name with screenshotNumber
  std::string screenshotNumberStr = std::to_string(screenshotNumber);
  int numZeros = 3 - screenshotNumberStr.size();
  std::string zeros(numZeros, '0');
  std::string filename = zeros + screenshotNumberStr + ".jpg";

  if (screenshotImg.save(filename.c_str(), ImageIO::FORMAT_JPEG) == ImageIO::OK) {
    cout << "File " << filename << " saved successfully." << endl;
  }
  else cout << "Failed to save file " << filename << '.' << endl;
}

void idleFunc()
{
  // animate coaster
  if (animationOn) {
    coaster.stepForward(animationSpeed);
  }
  // Notify GLUT that it should call displayFunc.
  glutPostRedisplay();
}

void reshapeFunc(int w, int h)
{
  glViewport(0, 0, w, h);

  // When the window has been resized, we need to re-set our projection matrix.
  matrix.SetMatrixMode(OpenGLMatrix::Projection);
  matrix.LoadIdentity();
  // You need to be careful about setting the zNear and zFar. 
  // Anything closer than zNear, or further than zFar, will be culled.
  const float zNear = 0.1f;
  const float zFar = 10000.0f;
  const float humanFieldOfView = 60.0f;
  matrix.Perspective(humanFieldOfView, 1.0f * w / h, zNear, zFar);
}

void mouseMotionDragFunc(int x, int y)
{
  // Mouse has moved, and one of the mouse buttons is pressed (dragging).

  // the change in mouse position since the last invocation of this function
  int mousePosDelta[2] = { x - mousePos[0], y - mousePos[1] };

  switch (controlState)
  {
    // translate the terrain
    case TRANSLATE:
      if (leftMouseButton)
      {
        // control x,y translation via the left mouse button
        terrainTranslate[0] += mousePosDelta[0] * 0.01f;
        terrainTranslate[1] -= mousePosDelta[1] * 0.01f;
      }
      if (middleMouseButton)
      {
        // control z translation via the middle mouse button
        terrainTranslate[2] += mousePosDelta[1] * 0.01f;
      }
      break;

    // rotate the terrain
    case ROTATE:
      if (leftMouseButton)
      {
        // control x,y rotation via the left mouse button
        terrainRotate[0] += mousePosDelta[1];
        terrainRotate[1] += mousePosDelta[0];
      }
      if (middleMouseButton)
      {
        // control z rotation via the middle mouse button
        terrainRotate[2] += mousePosDelta[1];
      }
      break;

    // scale the terrain
    case SCALE:
      if (leftMouseButton)
      {
        // control x,y scaling via the left mouse button
        terrainScale[0] *= 1.0f + mousePosDelta[0] * 0.01f;
        terrainScale[1] *= 1.0f - mousePosDelta[1] * 0.01f;
      }
      if (middleMouseButton)
      {
        // control z scaling via the middle mouse button
        terrainScale[2] *= 1.0f - mousePosDelta[1] * 0.01f;
      }
      break;
  }

  // store the new mouse position
  mousePos[0] = x;
  mousePos[1] = y;
}

void mouseMotionFunc(int x, int y)
{
  // Mouse has moved.
  // Store the new mouse position.
  mousePos[0] = x;
  mousePos[1] = y;
}

void mouseButtonFunc(int button, int state, int x, int y)
{
  // A mouse button has has been pressed or depressed.

  // Keep track of the mouse button state, in leftMouseButton, middleMouseButton, rightMouseButton variables.
  switch (button)
  {
    case GLUT_LEFT_BUTTON:
      leftMouseButton = (state == GLUT_DOWN);
    break;

    case GLUT_MIDDLE_BUTTON:
      middleMouseButton = (state == GLUT_DOWN);
    break;

    case GLUT_RIGHT_BUTTON:
      rightMouseButton = (state == GLUT_DOWN);
    break;
  }

  // Keep track of whether CTRL and SHIFT keys are pressed.
  switch (glutGetModifiers())
  {
    case GLUT_ACTIVE_CTRL:
      controlState = TRANSLATE;
    break;

    case GLUT_ACTIVE_SHIFT:
      controlState = SCALE;
    break;

    // If CTRL and SHIFT are not pressed, we are in rotate mode.
    default:
      controlState = ROTATE;
    break;
  }

  // Store the new mouse position.
  mousePos[0] = x;
  mousePos[1] = y;
}

void keyboardFunc(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 27: // ESC key
      exit(0); // exit the program
    break;

    case ' ':
      cout << "You pressed the spacebar." << endl;
    break;

    case 'x':
      // Take a screenshot.
      saveScreenshot(screenshotCount);
      ++screenshotCount;
    break;

    // speed up or slow down animation speed
    case '+':
    case '=':
      ++animationSpeed;
    break;
    case '-':
    case '_':
      --animationSpeed;
    break;

    case 's':
      animationOn = !animationOn;
    break;
  }
}

void displayFunc()
{
  // This function performs the actual rendering.

  // First, clear the screen.
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set up the camera position.
  matrix.SetMatrixMode(OpenGLMatrix::ModelView);
  matrix.LoadIdentity();

  // this method sets the camera at the current point,
  // looking along the tangent.
  coaster.lookAtOnCoaster();

  matrix.Translate(terrainTranslate[0], terrainTranslate[1], terrainTranslate[2]);
  matrix.Rotate(terrainRotate[0], 1.0, 0.0, 0.0);
  matrix.Rotate(terrainRotate[1], 0.0, 1.0, 0.0);
  matrix.Rotate(terrainRotate[2], 0.0, 0.0, 1.0);
  matrix.Scale(terrainScale[0], terrainScale[1], terrainScale[2]);
  // In here, you can do additional modeling on the object, such as performing translations, rotations and scales.
  // ...

  // Read the current modelview and projection matrices from our helper class.
  // The matrices are only read here; nothing is actually communicated to OpenGL yet.
  float modelViewMatrix[16];
  matrix.SetMatrixMode(OpenGLMatrix::ModelView);
  matrix.GetMatrix(modelViewMatrix);

  float projectionMatrix[16];
  matrix.SetMatrixMode(OpenGLMatrix::Projection);
  matrix.GetMatrix(projectionMatrix);

  // render groundWalls
  texturePipelineProgram.Bind();
  texturePipelineProgram.SetUniformVariableMatrix4fv("modelViewMatrix", GL_FALSE, modelViewMatrix);
  texturePipelineProgram.SetUniformVariableMatrix4fv("projectionMatrix", GL_FALSE, projectionMatrix);
  
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, groundWallsHandle);
  GLint h_groundWalls_textureImage = glGetUniformLocation(texturePipelineProgram.GetProgramHandle(), "textureImage");
  glUniform1i(h_groundWalls_textureImage, 0);
  
  vaogroundWalls.Bind();
  glDrawArrays(GL_TRIANGLES, 0, 30);

  // render sky
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, skyHandle);
  GLint h_sky_textureImage = glGetUniformLocation(texturePipelineProgram.GetProgramHandle(), "textureImage");
  glUniform1i(h_sky_textureImage, 0);
  
  vaoSky.Bind();
  glDrawArrays(GL_TRIANGLES, 0, 6);

  pipelineProgram.Bind();
  // Upload the modelview and projection matrices to the GPU. Note that these are "uniform" variables.
  // Important: these matrices must be uploaded to *all* pipeline programs used.
  // In hw1, there is only one pipeline program, but in hw2 there will be several of them.
  // In such a case, you must separately upload to *each* pipeline program.
  // Important: do not make a typo in the variable name below; otherwise, the program will malfunction.
  pipelineProgram.SetUniformVariableMatrix4fv("modelViewMatrix", GL_FALSE, modelViewMatrix);
  pipelineProgram.SetUniformVariableMatrix4fv("projectionMatrix", GL_FALSE, projectionMatrix);
  pipelineProgram.SetUniformVariableVec3f("lightPosition", lightPosition[0], lightPosition[1], lightPosition[2]);

  // Execute the rendering.
  // Bind the VAO that we want to render. Remember, one object = one VAO. 
  vao.Bind();
  eboRailIndeces.Bind();
  glDrawElements(GL_TRIANGLES, railIndices.size(), GL_UNSIGNED_INT, 0); // Render the VAO, by rendering "numVertices", starting from vertex 0.
  // Swap the double-buffers.
  glutSwapBuffers();
}

int initTexture(const char * imageFilename, GLuint textureHandle)
{
  // Read the texture image.
  ImageIO img;
  ImageIO::fileFormatType imgFormat;
  ImageIO::errorType err = img.load(imageFilename, &imgFormat);

  if (err != ImageIO::OK) 
  {
    printf("Loading texture from %s failed.\n", imageFilename);
    return -1;
  }

  // Check that the number of bytes is a multiple of 4.
  if (img.getWidth() * img.getBytesPerPixel() % 4) 
  {
    printf("Error (%s): The width*numChannels in the loaded image must be a multiple of 4.\n", imageFilename);
    cout << "Currently, it's " << img.getWidth() * img.getBytesPerPixel() << endl; 
    return -1;
  }

  // Allocate space for an array of pixels.
  int width = img.getWidth();
  int height = img.getHeight();
  unsigned char * pixelsRGBA = new unsigned char[4 * width * height]; // we will use 4 bytes per pixel, i.e., RGBA

  // Fill the pixelsRGBA array with the image pixels.
  memset(pixelsRGBA, 0, 4 * width * height); // set all bytes to 0
  for (int h = 0; h < height; h++)
    for (int w = 0; w < width; w++) 
    {
      // assign some default byte values (for the case where img.getBytesPerPixel() < 4)
      pixelsRGBA[4 * (h * width + w) + 0] = 0; // red
      pixelsRGBA[4 * (h * width + w) + 1] = 0; // green
      pixelsRGBA[4 * (h * width + w) + 2] = 0; // blue
      pixelsRGBA[4 * (h * width + w) + 3] = 255; // alpha channel; fully opaque

      // set the RGBA channels, based on the loaded image
      int numChannels = img.getBytesPerPixel();
      for (int c = 0; c < numChannels; c++) // only set as many channels as are available in the loaded image; the rest get the default value
        pixelsRGBA[4 * (h * width + w) + c] = img.getPixel(w, h, c);
    }

  // Bind the texture.
  glBindTexture(GL_TEXTURE_2D, textureHandle);

  // Initialize the texture.
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixelsRGBA);

  // Generate the mipmaps for this texture.
  glGenerateMipmap(GL_TEXTURE_2D);

  // Set the texture parameters.
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  // Query support for anisotropic texture filtering.
  GLfloat fLargest;
  glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &fLargest);
  printf("Max available anisotropic samples: %f\n", fLargest);
  // Set anisotropic texture filtering.
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 0.5f * fLargest);

  // Query for any errors.
  GLenum errCode = glGetError();
  if (errCode != 0) 
  {
    printf("Texture initialization error. Error code: %d.\n", errCode);
    return -1;
  }
  
  // De-allocate the pixel array -- it is no longer needed.
  delete [] pixelsRGBA;

  return 0;
}

void initTextures() {
  glGenTextures(1, &groundWallsHandle);
  glGenTextures(1, &skyHandle);
  int groundWallsInitStatus = initTexture(groundWallsImageFilePath, groundWallsHandle);
  int skyInitStatus = initTexture(skyImageFilePath, skyHandle);
  if (groundWallsInitStatus != 0 || skyInitStatus != 0) {
    cout << "Error initilizing textures in initTexture.\n";
    exit(-1);
  }

  const float groundWallsSkySize = 150.0f;
  const float groundWallsY = -5.0f;
  const float skyY = 30.0f;
  // ******************** init groundWalls: text coords, position, etc. ********************
  const float groundWallsVertices[] = {
    // ground
    -groundWallsSkySize, groundWallsY, -groundWallsSkySize,
    -groundWallsSkySize, groundWallsY, groundWallsSkySize,
    groundWallsSkySize, groundWallsY, -groundWallsSkySize,
    groundWallsSkySize, groundWallsY, groundWallsSkySize,
    -groundWallsSkySize, groundWallsY, groundWallsSkySize,
    groundWallsSkySize, groundWallsY, -groundWallsSkySize,

    // wall 1
    -groundWallsSkySize, groundWallsY, -groundWallsSkySize,
    groundWallsSkySize, groundWallsY, -groundWallsSkySize,
    groundWallsSkySize, skyY, -groundWallsSkySize,
    groundWallsSkySize, skyY, -groundWallsSkySize,
    -groundWallsSkySize, skyY, -groundWallsSkySize,
    -groundWallsSkySize, groundWallsY, -groundWallsSkySize,

    // wall 2
    -groundWallsSkySize, groundWallsY, groundWallsSkySize,
    -groundWallsSkySize, skyY, groundWallsSkySize,
    groundWallsSkySize, skyY, groundWallsSkySize,
    groundWallsSkySize, skyY, groundWallsSkySize,
    groundWallsSkySize, groundWallsY, groundWallsSkySize,
    -groundWallsSkySize, groundWallsY, groundWallsSkySize,

    // wall 3
    -groundWallsSkySize, groundWallsY, -groundWallsSkySize,
    -groundWallsSkySize, skyY, -groundWallsSkySize,
    -groundWallsSkySize, skyY, groundWallsSkySize,
    -groundWallsSkySize, skyY, groundWallsSkySize,
    -groundWallsSkySize, groundWallsY, groundWallsSkySize,
    -groundWallsSkySize, groundWallsY, -groundWallsSkySize,

    // wall 4
    groundWallsSkySize, groundWallsY, -groundWallsSkySize,
    groundWallsSkySize, groundWallsY, groundWallsSkySize,
    groundWallsSkySize, skyY, groundWallsSkySize,
    groundWallsSkySize, skyY, groundWallsSkySize,
    groundWallsSkySize, skyY, -groundWallsSkySize,
    groundWallsSkySize, groundWallsY, -groundWallsSkySize,
  };

  float groundWallsTexCoords[] = {
    // ground
    0.0f, 0.0f,
    1.0f, 0.0f,
    1.0f, 1.0f,
    1.0f, 1.0f,
    0.0f, 1.0f,
    0.0f, 0.0f,
    // wall 1
    0.0f, 0.0f,
    1.0f, 0.0f,
    1.0f, 1.0f,
    1.0f, 1.0f,
    0.0f, 1.0f,
    0.0f, 0.0f,
    // wall 2
    0.0f, 0.0f,
    1.0f, 0.0f,
    1.0f, 1.0f,
    1.0f, 1.0f,
    0.0f, 1.0f,
    0.0f, 0.0f,
    // wall 3
    0.0f, 0.0f,
    1.0f, 0.0f,
    1.0f, 1.0f,
    1.0f, 1.0f,
    0.0f, 1.0f,
    0.0f, 0.0f,
    // wall 4
    0.0f, 0.0f,
    1.0f, 0.0f,
    1.0f, 1.0f,
    1.0f, 1.0f,
    0.0f, 1.0f,
    0.0f, 0.0f,
  };

  vbogroundWallsVertices.Gen((6 * 5), 3, groundWallsVertices, GL_STATIC_DRAW);
  vbogroundWallsTexCoords.Gen((6*5), 2, groundWallsTexCoords, GL_STATIC_DRAW);
  vaogroundWalls.Gen();
  vaogroundWalls.ConnectPipelineProgramAndVBOAndShaderVariable(&texturePipelineProgram, &vbogroundWallsVertices, "position");
  vaogroundWalls.ConnectPipelineProgramAndVBOAndShaderVariable(&texturePipelineProgram, &vbogroundWallsTexCoords, "texCoord");
  
  // ******************** init sky: text coords, position, etc. ********************
  const float skyVertices[] = {
    -groundWallsSkySize, skyY, -groundWallsSkySize,
    -groundWallsSkySize, skyY, groundWallsSkySize,
    groundWallsSkySize, skyY, -groundWallsSkySize,
    groundWallsSkySize, skyY, groundWallsSkySize,
    -groundWallsSkySize, skyY, groundWallsSkySize,
    groundWallsSkySize, skyY, -groundWallsSkySize,
  };
  float skyTexCoords[] = {
    0.0f, 0.0f,
    1.0f, 0.0f,
    1.0f, 1.0f,
    1.0f, 1.0f,
    0.0f, 1.0f,
    0.0f, 0.0f
  };
  vboSkyVertices.Gen(6, 3, skyVertices, GL_STATIC_DRAW);
  vboSkyTexCoords.Gen(6, 2, skyTexCoords, GL_STATIC_DRAW);
  
  vaoSky.Gen();
  vaoSky.ConnectPipelineProgramAndVBOAndShaderVariable(&texturePipelineProgram, &vboSkyVertices, "position");
  vaoSky.ConnectPipelineProgramAndVBOAndShaderVariable(&texturePipelineProgram, &vboSkyTexCoords, "texCoord");
}

void initScene(int argc, char *argv[])
{
  // Set the backgroundWalls color.
  glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black color.

  // Enable z-buffering (i.e., hidden surface removal using the z-buffer algorithm).
  glEnable(GL_DEPTH_TEST);

  // Create a pipeline program. This operation must be performed BEFORE we initialize any VAOs.
  // A pipeline program contains our shaders. Different pipeline programs may contain different shaders.
  // In this homework, we only have one set of shaders, and therefore, there is only one pipeline program.
  // In hw2, we will need to shade different objects with different shaders, and therefore, we will have
  // several pipeline programs (e.g., one for the rails, one for the groundWalls/sky, etc.).
  // Load and set up the pipeline program, including its shaders.
  if (pipelineProgram.BuildShadersFromFiles(shaderBasePath, "vertexShader.glsl", "fragmentShader.glsl") != 0)
  {
    cout << "Failed to build the pipeline program." << endl;
    throw 1;
  } 
  cout << "Successfully built the pipeline program." << endl;
    
  // Bind the pipeline program that we just created. 
  // The purpose of binding a pipeline program is to activate the shaders that it contains, i.e.,
  // any object rendered from that point on, will use those shaders.
  // When the application starts, no pipeline program is bound, which means that rendering is not set up.
  // So, at some point (such as below), we need to bind a pipeline program.
  // From that point on, exactly one pipeline program is bound at any moment of time.
  pipelineProgram.Bind();

  // Prepare the triangle position and color data for the VBO. 
  // The code below sets up a single triangle (3 vertices).
  // The triangle will be rendered using GL_TRIANGLES (in displayFunc()).

  // Create the VBOs. 
  // We make a separate VBO for vertices and normals. 
  // This operation must be performed BEFORE we initialize any VAOs.
  vboVertices.Gen(numRailPoints, 3, railPoints.data(), GL_STATIC_DRAW); // 3 values per position
  eboRailIndeces.Gen(railIndices.size() * sizeof(unsigned int), (const int*) railIndices.data(), GL_STATIC_DRAW);
  vboNormals.Gen(numRailPoints, 3, railNormals.data(), GL_STATIC_DRAW); // 3 values per normal

  // Create the VAOs. There is a single VAO in this example.
  // Important: this code must be executed AFTER we created our pipeline program, and AFTER we set up our VBOs.
  // A VAO contains the geometry for a single object. There should be one VAO per object.
  // In this homework, "geometry" means vertex positions and normals. In homework 2, it will also include
  // vertex normal and vertex texture coordinates for texture mapping.
  vao.Gen();

  // Set up the relationship between the "position" shader variable and the VAO.
  // Important: any typo in the shader variable name will lead to malfunction.
  vao.ConnectPipelineProgramAndVBOAndShaderVariable(&pipelineProgram, &vboVertices, "position");

  // Set up the relationship between the "normal" shader variable and the VAO.
  // Important: any typo in the shader variable name will lead to malfunction.
  vao.ConnectPipelineProgramAndVBOAndShaderVariable(&pipelineProgram, &vboNormals, "normal");

  // texture stuff:
  if (texturePipelineProgram.BuildShadersFromFiles(shaderBasePath, "textureVertexShader.glsl", "textureFragmentShader.glsl") != 0)
  {
    cout << "Failed to build the texture pipeline program." << endl;
    throw 1;
  }
  cout << "Successfully built the texture pipeline program." << endl;
  texturePipelineProgram.Bind();
  initTextures();

  // Check for any OpenGL errors.
  std::cout << "GL error status is: " << glGetError() << std::endl;
}

// Given a fully populated vector of points along a spline,
// populates "points" data (for VBO) and "indeces" data (for EBO)
// to render the coaster track. 
void initRail() {
  unsigned int numRailChunks = numVertices - 1;
  unsigned int numPointsPerChunk = 8;
  railPoints.resize(numRailChunks * numPointsPerChunk * 3, 0.0f);
  numRailPoints = numRailChunks * numPointsPerChunk;
  railNormals.resize(numRailPoints * 3, 0.0f);
  unsigned int railIdx = 0;
  unsigned int base = 0;
  
  for (int i = 0; i < points.size() - 6; i += 3) {
    // get points
    float p1[3] = {points[i], points[i+1], points[i+2]};
    float p2[3] = {points[i+3], points[i+4], points[i+5]};
    float t1[3] = {tangents[i], tangents[i+1], tangents[i+2]};
    float t2[3] = {tangents[i+3], tangents[i+4], tangents[i+5]};
    
    // find normals
    float n1[3], n2[3];
    cross(t1, worldUp, n1);
    normalize(n1);
    cross(t2, worldUp, n2);
    normalize(n2);
    
    // find binormals
    float b1[3], b2[3];
    cross(t1, n1, b1);
    normalize(b1);
    cross(t2, n2, b2);
    normalize(b2);
    
    // compute vertex positions
    for (int j = 0; j < 3; ++j) {
        railPoints[railIdx + j] = p1[j] + alpha * (-b1[j] - n1[j]);
        railPoints[railIdx + j + 3] = p1[j] + alpha * (b1[j] - n1[j]);
        railPoints[railIdx + j + 6] = p1[j] + alpha * (b1[j] + n1[j]);
        railPoints[railIdx + j + 9] = p1[j] + alpha * (-b1[j] + n1[j]);
        railPoints[railIdx + j + 12] = p2[j] + alpha * (-b2[j] - n2[j]);
        railPoints[railIdx + j + 15] = p2[j] + alpha * (b2[j] - n2[j]);
        railPoints[railIdx + j + 18] = p2[j] + alpha * (b2[j] + n2[j]);
        railPoints[railIdx + j + 21] = p2[j] + alpha * (-b2[j] + n2[j]);
    }
    
    // find per vertex normals
    // v0
    float norm0[3] = {-b1[0] - n1[0], -b1[1] - n1[1], -b1[2] - n1[2]};
    normalize(norm0);
    railNormals[railIdx + 0] = norm0[0];
    railNormals[railIdx + 1] = norm0[1];
    railNormals[railIdx + 2] = norm0[2];

    // v1
    float norm1[3] = {b1[0] - n1[0], b1[1] - n1[1], b1[2] - n1[2]};
    normalize(norm1);
    railNormals[railIdx + 3] = norm1[0];
    railNormals[railIdx + 4] = norm1[1];
    railNormals[railIdx + 5] = norm1[2];
    
    // v2
    float norm2[3] = {b1[0] + n1[0], b1[1] + n1[1], b1[2] + n1[2]};
    normalize(norm2);
    railNormals[railIdx + 6] = norm2[0];
    railNormals[railIdx + 7] = norm2[1];
    railNormals[railIdx + 8] = norm2[2];

    // v3
    float norm3[3] = {-b1[0] + n1[0], -b1[1] + n1[1], -b1[2] + n1[2]};
    normalize(norm3);
    railNormals[railIdx + 9] = norm3[0];
    railNormals[railIdx + 10] = norm3[1];
    railNormals[railIdx + 11] = norm3[2];
    
    // v4
    float norm4[3] = {-b2[0] - n2[0], -b2[1] - n2[1], -b2[2] - n2[2]};
    normalize(norm4);
    railNormals[railIdx + 12] = norm4[0];
    railNormals[railIdx + 13] = norm4[1];
    railNormals[railIdx + 14] = norm4[2];

    // v5
    float norm5[3] = {b2[0] - n2[0], b2[1] - n2[1], b2[2] - n2[2]};
    normalize(norm5);
    railNormals[railIdx + 15] = norm5[0];
    railNormals[railIdx + 16] = norm5[1];
    railNormals[railIdx + 17] = norm5[2];
    
    // v6
    float norm6[3] = {b2[0] + n2[0], b2[1] + n2[1], b2[2] + n2[2]};
    normalize(norm6);
    railNormals[railIdx + 18] = norm6[0];
    railNormals[railIdx + 19] = norm6[1];
    railNormals[railIdx + 20] = norm6[2];

    // v7
    float norm7[3] = {-b2[0] + n2[0], -b2[1] + n2[1], -b2[2] + n2[2]};
    normalize(norm7);
    railNormals[railIdx + 21] = norm7[0];
    railNormals[railIdx + 22] = norm7[1];
    railNormals[railIdx + 23] = norm7[2];
    
    // top square
    drawSquare(base + 1, base + 2, base + 5, base + 6, railIndices);
    // left square
    drawSquare(base + 2, base + 3, base + 6, base + 7, railIndices);
    // bottom square
    drawSquare(base + 3, base + 0, base + 7, base + 4, railIndices);
    // right square
    drawSquare(base + 0, base + 1, base + 4, base + 5, railIndices);
    
    railIdx += numPointsPerChunk * 3;
    base += numPointsPerChunk;
  }
}

void initCoaster(char * argv) 
{
  FILE * fileSpline = fopen(argv, "r");
  if (fileSpline == NULL) 
  {
    printf ("Cannot open file %s.\n", argv);
    exit(1);
  }

  // Read the number of spline control points.
  fscanf(fileSpline, "%d\n", &spline.numControlPoints);
  printf("Detected %d control points.\n", spline.numControlPoints);

  // allocate memory.
  spline.points.resize(spline.numControlPoints);
  // Load the control points.
  for(int i=0; i<spline.numControlPoints; i++)
  {
    if (fscanf(fileSpline, "%f %f %f", 
           &spline.points[i].x, 
	   &spline.points[i].y, 
	   &spline.points[i].z) != 3)
    {
      printf("Error: incorrect number of control points in file %s.\n", argv);
      exit(1);
    }
  }

  // ******************* calculate points along spline *******************
  float parameterVector[4];
  vector<vector<float>> blendingWeights;

  for (float u = 0.00; u < 1.0; u += 0.001)
  {
    // param vector [u^3, u^2, u, 1]
    parameterVector[0] = u * u * u;
    parameterVector[1] = u * u;
    parameterVector[2] = u;
    parameterVector[3] = 1.0f;

    // compute blended weights = basisMatrix * parameterVector
    float blendedWeight[4];
    MultiplyMatrices(4, 4, 1, basisMatrix, parameterVector, blendedWeight);

    blendingWeights.push_back({ blendedWeight[0], blendedWeight[1],
      blendedWeight[2], blendedWeight[3] });
  }

  // control matrix per segment.
  // 4 x 4
  // (p1, p2, p3, p4)
  vector<vector<float>> pointControlMatrix;

  for (int i = 1; i < spline.numControlPoints - 2; ++i)
  {
    pointControlMatrix = {
      {spline.points[i-1].x, spline.points[i-1].y, spline.points[i-1].z},
      {spline.points[i].x, spline.points[i].y, spline.points[i].z},
      {spline.points[i+1].x, spline.points[i+1].y, spline.points[i+1].z},
      {spline.points[i+2].x, spline.points[i+2].y, spline.points[i+2].z}
    };
    for (auto &bw : blendingWeights)
    {
      // x
      points.push_back(
      bw[0] * pointControlMatrix[0][0] +
      bw[1] * pointControlMatrix[1][0] +
      bw[2] * pointControlMatrix[2][0] +
      bw[3] * pointControlMatrix[3][0]);

      // y
      points.push_back(
      bw[0] * pointControlMatrix[0][1] +
      bw[1] * pointControlMatrix[1][1] +
      bw[2] * pointControlMatrix[2][1] +
      bw[3] * pointControlMatrix[3][1]);

      // z
      points.push_back(
      bw[0] * pointControlMatrix[0][2] +
      bw[1] * pointControlMatrix[1][2] +
      bw[2] * pointControlMatrix[2][2] +
      bw[3] * pointControlMatrix[3][2]);
    }
  }
  fclose(fileSpline);
  printf("Loaded spline with %d control point(s).\n", spline.numControlPoints);
  numVertices = points.size()/3; // This must be a global variable, so that we know how many vertices to render in glDrawArrays.
  // ****** compute tangents for every point along interpolated spline. ******
  float u = 0.00;
  float tagnent[3];
  // most of the tangents done here.
  for (int i = 0; i < points.size(); i += 3) {
    float uMatrix[4] = {3 * (u*u), 2*u, 1, 0};
    float temp[4];
    float tangentControlMatrix[12];
    // dynamically figure out control matrix
    for (int j = -1; j <= 2; ++j) {
      int idx = i + j * 3;

      // clamp so we dont go out of bounds
      if (idx < 0) idx = 0;
      if (idx > static_cast<int>(points.size() - 3)) {
        idx = points.size() - 3;
      }

      tangentControlMatrix[j + 1] = points[idx];
      tangentControlMatrix[j + 5] = points[idx + 1];
      tangentControlMatrix[j + 9] = points[idx + 2];
    }
    MultiplyMatrices(4, 4, 1, basisMatrix, uMatrix, temp);
    // x
    tagnent[0] = temp[0] * tangentControlMatrix[0] + temp[1] * tangentControlMatrix[1] + 
      temp[2] * tangentControlMatrix[2] + temp[3] * tangentControlMatrix[3];
    // y
    tagnent[1] = temp[0] * tangentControlMatrix[4] + temp[1] * tangentControlMatrix[5] + 
      temp[2] * tangentControlMatrix[6] + temp[3] * tangentControlMatrix[7];
    // z
    tagnent[2] = temp[0] * tangentControlMatrix[8] + temp[1] * tangentControlMatrix[9] + 
      temp[2] * tangentControlMatrix[10] + temp[3] * tangentControlMatrix[11];
    normalize(tagnent);
    if (u >= 1.0f) u = 0.0f;
    else u+=0.001f;
    tangents.push_back(tagnent[0]);
    tangents.push_back(tagnent[1]);
    tangents.push_back(tagnent[2]);
  }
  // ******************* spline done. next, generate the rail. ***************
  initRail();
}

void initOpenGl(int argc, char *argv[]) {

  cout << "Initializing GLUT..." << endl;
  glutInit(&argc,argv);

  cout << "Initializing OpenGL..." << endl;

  #ifdef __APPLE__
    glutInitDisplayMode(GLUT_3_2_CORE_PROFILE | GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
  #else
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
  #endif

  glutInitWindowSize(windowWidth, windowHeight);
  glutInitWindowPosition(0, 0);  
  glutCreateWindow(windowTitle);

  cout << "OpenGL Version: " << glGetString(GL_VERSION) << endl;
  cout << "OpenGL Renderer: " << glGetString(GL_RENDERER) << endl;
  cout << "Shading Language Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;

  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(windowWidth - 1, windowHeight - 1);
  #endif

  // Tells GLUT to use a particular display function to redraw.
  glutDisplayFunc(displayFunc);
  // Perform animation inside idleFunc.
  glutIdleFunc(idleFunc);
  // callback for mouse drags
  glutMotionFunc(mouseMotionDragFunc);
  // callback for idle mouse movement
  glutPassiveMotionFunc(mouseMotionFunc);
  // callback for mouse button changes
  glutMouseFunc(mouseButtonFunc);
  // callback for resizing the window
  glutReshapeFunc(reshapeFunc);
  // callback for pressing the keys on the keyboard
  glutKeyboardFunc(keyboardFunc);

  // init glew
  #ifdef __APPLE__
    // nothing is needed on Apple
  #else
    // Windows, Linux
    GLint result = glewInit();
    if (result != GLEW_OK)
    {
      cout << "error: " << glewGetErrorString(result) << endl;
      exit(EXIT_FAILURE);
    }
  #endif
}

// Note: You should combine this file with the solution of homework 1.

// Note for Windows/MS Visual Studio:
// You should set argv[1] to the name of your spline, for example, splines/circle.sp .
// To do this, on the "Solution Explorer", right click your project, choose "Properties",
// go to "Configuration Properties", click "Debug",
// then type your track file name for the "Command Arguments".
// You can also repeat this process for the "Release" configuration.

int main (int argc, char ** argv)
{
  if (argc < 2)
  {  
    printf ("Usage: %s <spline file>\n", argv[0]);
    exit(0);
  }

  // This function loads the spline from the provided filename,
  // interpolates the spline, and then generates the points for the
  // rail based on these interpolated points.
  initCoaster(argv[1]);

  // init opengl
  initOpenGl(argc, argv);

  // Perform the initialization.
  // Builds all shader programs.
  // Creates VBOs and VAOs.
  // note: the latter part of this function calls
  // initTextures, which renders the geometry
  // for the skybox/ground/walls and creates the
  // UV maps needed to use the textures.
  initScene(argc, argv);

  // Sink forever into the GLUT loop.
  glutMainLoop();

  return 0;
}