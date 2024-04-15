#include "GLout.h"

void processInput(GLFWwindow *window)
{
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

void GLout(Mat img,vector<vector<Point> > mesh,vector<vector<Point> > V)
{
    int n = img.rows,m = img.cols;
    
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); 

    GLFWwindow*window = glfwCreateWindow(m,n,"Out Resizing",NULL,NULL); //创建一个窗口对象

    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) //初始化GLAD
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        exit(0);
    }

    Shader ourShader("E:/c++/wl/src/texture/texture.vs", "E:/c++/wl/src/texture/texture.fx"); 
    
    // float vertices[] = {
    //     // positions          // colors           // texture coords
    //      0.5f,  0.5f, 0.0f,   0.0f, 0.0f, 0.0f,   1.0f, 1.0f, // top right
    //      0.5f, -0.5f, 0.0f,   0.0f, 0.0f, 0.0f,   1.0f, 0.0f, // bottom right
    //     -0.5f, -0.5f, 0.0f,   0.0f, 0.0f, 0.0f,   0.0f, 0.0f, // bottom left
    //     -0.5f,  0.5f, 0.0f,   0.0f, 0.0f, 0.0f,   0.0f, 1.0f, // top left 
    //      0.0f,  0.0f, 0.0f,   0.0f, 0.0f, 0.0f,   0.5f, 0.5f  // center
    // };
    // unsigned int indices[] = {  
    //     0, 1, 3, // first triangle
    //     0, 3, 4  // second triangle
    // };
    // float vertices[] = {
    //     // positions          // colors           // texture coords
    //      -1.0f, -1.0f, 0.0f,   0.0f, 0.0f, 0.0f,   0.02698413f, 0.11498258f, // bottom left
    //       0.99365079,-1.00000000, 0.0f,   0.0f, 0.0f, 0.0f,   0.94126981, 0.00000000,  // bottom right
    //       0.99365079,0.98606277, 0.0f,   0.0f, 0.0f, 0.0f,   0.99206346,0.81533098  // bottom right
    // };
    // unsigned int indices[] = {  
    //     0, 1, 21, // first triangle
    // };
    // unsigned int indices[] = {1,20,20*21+20};
    
    float vertices[21 * 21 * 8];
    // float*vertices = new float[21 * 21 * 8];
    for(int i = 0;i <= 20; ++ i)
        for(int j = 0;j <= 20; ++ j)
        {
            int ii = 20 - i,jj = j;
            // int x = (n - 1) / 20.0 * ii;
            // if(x > n - 1) x = n - 1;
            // int y = (m - 1) / 20.0 * jj;
            // if(y > m - 1) y = m - 1;
            int x = V[ii][jj].x,y = V[ii][jj].y;
            // vertices[(ii*21+jj)*8] = (float)2.0*mesh[20-ii][jj].y/m - 1,
            // vertices[(ii*21+jj)*8 + 1] = (float)2.0*mesh[20-ii][jj].x/n - 1;
            vertices[(i*21+j)*8] = (float)2.0*y/m - 1,
            vertices[(i*21+j)*8 + 1] = (float)2.0*x/n - 1;
            vertices[(i*21+j)*8 + 2] = (float)0;
            
            // vertices[(i*21+j)*8 + 3] = (float)img.at<Vec3b>(mesh[ii][jj].x,mesh[ii][jj].y)[2]/255.0,
            // vertices[(i*21+j)*8 + 4] = (float)img.at<Vec3b>(mesh[ii][jj].x,mesh[ii][jj].y)[1]/255.0;
            // vertices[(i*21+j)*8 + 5] = (float)img.at<Vec3b>(mesh[ii][jj].x,mesh[ii][jj].y)[0]/255.0;
            
            // vertices[(ii*21+jj)*8 + 6] = (float)1.0*y/m;
            // vertices[(ii*21+jj)*8 + 7] = (float)1.0*x/n;
            vertices[(i*21+j)*8 + 6] = (float)1.0*mesh[i][j].y/m-1;
            vertices[(i*21+j)*8 + 7] = (float)1.0*mesh[i][j].x/n-1;
        }

    // unsigned int*indices = new unsigned int[20*20*6];
    unsigned int indices[20*20*6];
    for(int i = 0;i < 20; ++ i)
        for(int j = 0;j < 20; ++ j)
        {
            indices[(i*20+j)*6] = i*21+j;
            indices[(i*20+j)*6+1] = i*21+j+1;
            indices[(i*20+j)*6+2] = (i+1)*21+j;
            indices[(i*20+j)*6+3] = i*21+j+1;
            indices[(i*20+j)*6+4] = (i+1)*21+(j+1);
            indices[(i*20+j)*6+5] = (i+1)*21+j;
        }
    // for(int i = 0;i < 3; ++ i)
    // {
    //     cout <<"indices: ";
    //     cout << indices[i] << ' ';cout << endl;
    //     cout << "vertices: ";
    //     for(int j = 0;j < 8;++ j)
    //     {
    //         cout <<fixed << setprecision(8) << vertices[indices[i]*8 + j] << ' ';
    //     }
    //     cout << endl;
    // }

    unsigned int VAO,VBO,EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    
    // 1. 绑定VAO,VBO
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);  
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
    // 2. 把顶点数组复制到缓冲中供OpenGL使用
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    // 3. 设置顶点属性指针
    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    // texture coord attribute
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);

    GLuint texture;
    glGenTextures(1, &texture);
    
    // 绑定这个纹理
    glBindTexture(GL_TEXTURE_2D, texture);

    // 设置纹理的缩放和环绕方式
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // OpenCV加载的图像通常是BGR格式，OpenGL需要的是RGB格式
    // cout << "buffer solve" << endl;
    // unsigned char* buffer = new unsigned char[n*m*3];
    // for(int i = 0;i < n; ++ i)
    //     for(int j = 0;j < m;++ j)
    //     {
    //         // cout <<(j*n + i) * 3 + 2 << " " << n * m * 3 << endl;
    //         buffer[(j + i*m) * 3] = img.at<Vec3b>(i,j)[0],
    //         buffer[(j + i*m) * 3 + 1] = img.at<Vec3b>(i,j)[1],
    //         buffer[(j + i*m) * 3 + 2] = img.at<Vec3b>(i,j)[2];
    //     }
    // cout << "buffer solve out" << endl;
    // // // 将图像数据上传到纹理
    Mat tempMat;
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); //宽为奇数
    cvtColor(img, tempMat, COLOR_BGR2RGB);
    flip(tempMat, tempMat, 0);
    flip(tempMat, tempMat, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, m, n, 0, GL_RGB, GL_UNSIGNED_BYTE, tempMat.data);

    // int width, height, nrChannels;
    // stbi_set_flip_vertically_on_load(true);
    // unsigned char *data = stbi_load("E:/c++/resizing/output.bmp",&width, &height, &nrChannels, 0);
    // unsigned char *data = stbi_load("E:/c++/resizing/input/3_input.jpg",&width, &height, &nrChannels, 0);
    // // cout << width << ' ' << height << endl;
    // if (data)
    // {
    //     glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);
    // }
    // else
    // {
    //     std::cout << "Failed to load texture" << std::endl;
    // }
    // stbi_image_free(data);
    // 生成 mipmaps，如果需要的话
    // glGenerateMipmap(GL_TEXTURE_2D);

    while(!glfwWindowShouldClose(window))
    {

        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // bind Texture
        glBindTexture(GL_TEXTURE_2D, texture);

        // render container
        ourShader.use();
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, 6*400 , GL_UNSIGNED_INT, 0);

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);

    glfwTerminate();
    return;
}