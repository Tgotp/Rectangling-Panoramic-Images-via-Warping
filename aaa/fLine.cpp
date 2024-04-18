#include "fLine.h"
double len_vector(point a) { return sqrt(a.x*a.x + a.y * a.y); }
double* fline::line_img(Mat img,int *n_line,Mat&line_img,bool is_show)
{
    int n = img.rows,m = img.cols;
    double *input_img = solve_img(img);
    double *out = lsd(n_line,input_img,m,n);
    
    img.copyTo(line_img);
    // line_img = ~line_img;
    int num = * n_line;
    for(int i = 0;i < num;++ i)
    {
        // if(len_vector(point(out[i*7],out[i*7 + 1]) - point(out[i*7 + 2],out[i*7 + 3])) < 100) continue;
        // cout << "Limited line"<< out[i*7]<< ' ' << out[i*7+1] << ' ' << out[i*7+2] << ' ' << out[i*7+3] << endl;
        line(line_img,Point(out[i*7],out[i*7+1]),Point(out[i*7+2],out[i*7+3]),Scalar(0,0,255),1);
    }
    if(is_show)
    {
        namedWindow("line",WINDOW_FREERATIO);
        imshow("line",line_img);
        waitKey(0);
    }
    return out;
}

bool onSegment(point p, point q, point r) {
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
        return true;
    return false;
}

int orientation(point p, point q, point r) {
    double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (val == 0) return 0; 
    return (val > 0) ? 1 : 2; 
}

bool doIntersect(Line l1, Line l2) {
    int o1 = orientation(l1.x, l1.y, l2.x);
    int o2 = orientation(l1.x, l1.y, l2.y);
    int o3 = orientation(l2.x, l2.y, l1.x);
    int o4 = orientation(l2.x, l2.y, l1.y);

    if (o1 != o2 && o3 != o4)
        return true;

    if (o1 == 0 && onSegment(l1.x, l2.x, l1.y)) return true;
    if (o2 == 0 && onSegment(l1.x, l2.y, l1.y)) return true;
    if (o3 == 0 && onSegment(l2.x, l1.x, l2.y)) return true;
    if (o4 == 0 && onSegment(l2.x, l1.y, l2.y)) return true;

    return false;
}
point lineIntersection(Line l1, Line l2) {
    // a1x + b1y = c1
    double a1 = l1.y.y - l1.x.y;
    double b1 = l1.x.x - l1.y.x;
    double c1 = a1*(l1.x.x) + b1*(l1.x.y);

    // a2x + b2y = c2
    double a2 = l2.y.y - l2.x.y;
    double b2 = l2.x.x - l2.y.x;
    double c2 = a2*(l2.x.x) + b2*(l2.x.y);

    double determinant = a1*b2 - a2*b1;

    int x = (b2*c1 - b1*c2) / determinant;
    int y = (a1*c2 - a2*c1) / determinant;
    return point(x, y);
}

Mat lineImg;
Line getLine(Line x,vector<Point> p,int z,bool&k)
{
    point t[5];int cnt = 0;k = z;
    for(int i = 0;i <= 3;++ i)
    {
        Line y = Line(point(p[i].x,p[i].y),point(p[(i+1)%4].x,p[(i+1)%4].y));
        
        if(doIntersect(x, y))
        {
            point L = lineIntersection(x, y);
            if(z)
            {
                if(L.x == INT_MIN && L.y == INT_MIN)
                    return x;
                return z == 1 ? Line(L,x.x) : Line(L,x.y);
            }
            else t[++ cnt] = L;
        }
    }
    k = cnt >= 2;
    // cout << k << endl;
    if(!k) return x;
    else return Line(t[1],t[2]);
}

vector<vector<vector<Line > > > fline::init_line(vector<vector<point > > mesh,double *seg,int n,int *line_n,int length)
{
    // cout << "init line" << endl;
    vector<vector<vector<Line > > > Mesh_Line;
    vector<vector<Line > > Mesh_line;
    vector<Line > mesh_line;
    vector<Point > contours;
    *line_n = 0;
    
    for(int i = 0;i < 20; ++ i)
    {
        Mesh_line.clear();
        for(int j = 0;j < 20;++ j)
        {
            mesh_line.clear();
            Point lt = Point(mesh[i][j].y,mesh[i][j].x);
            Point rt = Point(mesh[i][j+1].y,mesh[i][j+1].x);
            Point rb = Point(mesh[i+1][j+1].y,mesh[i+1][j+1].x);
            Point lb = Point(mesh[i+1][j].y,mesh[i+1][j].x);
            contours.clear();
            contours.push_back(lt);
            contours.push_back(rt);
            contours.push_back(rb);
            contours.push_back(lb);
            for(int k = 0;k < n;++ k) 
            {
                point st = point(seg[k*7],seg[k*7+1]);
                point ed = point(seg[k*7+2],seg[k*7+3]);
                if(len_vector(st - ed) < length) continue;
                int z1 = pointPolygonTest(contours,Point2f(st.x,st.y),false);
                int z2 = pointPolygonTest(contours,Point2f(ed.x,ed.y),false);
                bool t;
                if(z1 >= 0 && z2 >= 0)
                {
                    if(st == ed) continue;
                    ++ (*line_n);
                    mesh_line.push_back(Line(st,ed));
                }
                else
                {
                    Line a = getLine(Line(st,ed),contours,z1>=0 ? 1 : z2 >= 0? 2 : 0,t);
                    if(a.x == a.y || !t) continue;
                    ++ (*line_n);
                    mesh_line.push_back(a);
                }
            }
            Mesh_line.push_back(mesh_line);
        }
        Mesh_Line.push_back(Mesh_line);
    }

    // cout << *line_n << endl;
    // cout << "init line end" << endl;
    return Mesh_Line;
}
void fline::init_bins(int *num,vector<vector<vector<Line > > >&mesh_line,int number)
{
    for(int i = 0;i < number; ++ i) num[i] = 0;
    for(int i = 0;i < 20;++ i)
        for(int j = 0;j < 20;++ j)
            for(auto&k: mesh_line[i][j])
            {
                point q = point(k.x.x - k.y.x,k.x.y - k.y.y);
                k.theta = q.x / (sqrt(q.x * q.x + q.y * q.y));
                if(k.theta > 1) // too big
                    k.theta = 1;
                k.theta = acos(k.theta);
                k.pos = floor(k.theta/PI * number);
                k.pos = min(number - 1,k.pos);
                num[k.pos] ++ ;
                // cout << k.theta << endl;
            }
}
void fline::check_line(Mat img,vector<vector<vector<Line > > >mesh_line)
{
    cout << "check line" << endl;
    Mat fimg;
    img.copyTo(fimg);
    for(int i = 0;i < 20;++ i)
        for(int j = 0;j < 20;++ j)
            for(auto k : mesh_line[i][j])
            {
                
                // line(fimg,k.x,k.y,Scalar(0,255,0),1);
                line(fimg,Point(k.x.x,k.x.y),Point(k.y.x,k.y.y),Scalar(rand()%255,rand()%255,rand()%255),1);
                // namedWindow("mesh line",WINDOW_FREERATIO);
                // imshow("mesh line",fimg);
                // waitKey(0);
            }
        
    namedWindow("mesh line",WINDOW_FREERATIO);
    imshow("mesh line",fimg);
    waitKey(0);
}

double fline::get_gray(Mat&img,int i,int j) { return img.at<Vec3b>(i,j)[0] * 0.299 + img.at<Vec3b>(i,j)[1] * 0.587 + img.at<Vec3b>(i,j)[2] * 0.114; }

double* fline::solve_img(Mat img)
{
    int n = img.rows,m = img.cols;
    double*ans = new double[n * m];
    for(int j = 0;j < m; ++ j)
        for(int i = 0;i < n; ++ i)
        {
            // cout << n << ' ' << m <<' '<< i << ' ' << j <<' ' <<i +j*n <<endl;
            ans[i*m+j] = get_gray(img,i,j);
        }
    return ans;
}
