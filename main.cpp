/*****************************************************************************************
Author: João Fontes Gonçalves

Implementation of an heuristic for Minimum/Maximum Area Polygonizations

To compile: g++ -o main main.cpp -I/usr/include/python2.7 -lpython2.7

*****************************************************************************************/

#include <bits/stdc++.h>
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

const double eps = 1e-8;
const double PI = acos(-1.0);

struct pt {
    double x, y;
};

///////////////////////////////////////////////////

struct Triangle{
	pt a, b, c;
};

///////////////////////////////////////////////////

bool Equal(double a, double b){
	return abs(a-b) < eps;
}

///////////////////////////////////////////////////

bool Equal_points(pt a, pt b){
	return Equal(a.x, b.x) && Equal(a.y, b.y);
}

///////////////////////////////////////////////////

bool cmp(pt a, pt b) {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
}

///////////////////////////////////////////////////

bool cw(pt a, pt b, pt c) {
    return a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y) < 0;
}

///////////////////////////////////////////////////

bool ccw(pt a, pt b, pt c) {
    return a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y) > 0;
}

///////////////////////////////////////////////////

double dist(pt a, pt b){
	return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

///////////////////////////////////////////////////

double norm(pt a){
	return sqrt(a.x*a.x + a.y*a.y);
}

///////////////////////////////////////////////////

double dot_product(pt o, pt a, pt b){
	pt oa, ob;
	oa.x = a.x - o.x;
	oa.y = a.y - o.y;
	ob.x = b.x - o.x;
	ob.y = b.y - o.y;
	return oa.x*ob.x + oa.y*ob.y;
}

///////////////////////////////////////////////////

double cos_angle(pt o, pt a, pt b){
	pt oa, ob;
	oa.x = a.x - o.x;
	oa.y = a.y - o.y;
	ob.x = b.x - o.x;
	ob.y = b.y - o.y;
	return dot_product(o, a, b)/(norm(oa)*norm(ob));
}

///////////////////////////////////////////////////

//Test if points p and q are on the same side of line ab
bool same_side(pt a, pt b, pt p, pt q){
	return (cw(a,p,b) && cw(a,q,b)) || (ccw(a,p,b) && ccw(a,q,b));
}

///////////////////////////////////////////////////

//Test if point p is on segment ab
bool point_in_segment(pt p, pt a, pt b){
	return Equal(dist(p,a) + dist(p,b), dist(a,b));
}

///////////////////////////////////////////////////

//Test if point p is on line ab
bool point_in_line(pt p, pt a, pt b){
	return Equal(dist(p,a) + dist(p,b), dist(a,b)) || Equal(abs(dist(p,a) - dist(p,b)), dist(a,b));
}

///////////////////////////////////////////////////

//Test for intersection between two lines ab and cd
//Supposing that the two lines are never coincidents
bool line_intersect(pt a, pt b, pt c, pt d, pt &intersect){
	double m1, h1, m2, h2;

	if(Equal(a.x, b.x) && Equal(c.x, d.x)) return false;
	else if(Equal(a.x, b.x) && !Equal(c.x, d.x)){
		m2 = (d.y - c.y)/(d.x - c.x);
		h2 = d.y - m2*d.x;
		intersect.x = a.x;
		intersect.y = m2*a.x + h2;
		return true;
	}
	else if(!Equal(a.x, b.x) && Equal(c.x, d.x)){
		m1 = (b.y - a.y)/(b.x - a.x);
		h1 = b.y - m1*b.x;
		intersect.x = c.x;
		intersect.y = m1*c.x + h1;
		return true;
	}
	else{
		m1 = (b.y - a.y)/(b.x - a.x);
		h1 = b.y - m1*b.x;
		m2 = (d.y - c.y)/(d.x - c.x);
		h2 = d.y - m2*d.x;

		if(Equal(m1, m2)) return false;
		else{
			intersect.x = (h2 - h1)/(m1 - m2);
			intersect.y = m1*intersect.x + h1;
			return true;
		}
	}
}

///////////////////////////////////////////////////

//Test for intersection between line ab and segment cd
bool line_intersect_segment(pt a, pt b, pt c, pt d, pt &intersect){
	return line_intersect(a, b, c, d, intersect) && point_in_segment(intersect, c, d);
}

///////////////////////////////////////////////////

//Test for intersection between segment ab and segment cd
bool segment_intersect_segment(pt a, pt b, pt c, pt d, pt &intersect){
	return line_intersect(a, b, c, d, intersect) && point_in_segment(intersect, a, b) && point_in_segment(intersect, c, d);
}

///////////////////////////////////////////////////

//Get convex hull in clock-wise sense
void convex_hull(vector<pt>& a) {
    if (a.size() == 1)
        return;

    sort(a.begin(), a.end(), &cmp);
    pt p1 = a[0], p2 = a.back();
    vector<pt> up, down;
    up.push_back(p1);
    down.push_back(p1);
    for (int i = 1; i < (int)a.size(); i++) {
        if (i == a.size() - 1 || cw(p1, a[i], p2)) {
            while (up.size() >= 2 && !cw(up[up.size()-2], up[up.size()-1], a[i]))
                up.pop_back();
            up.push_back(a[i]);
        }
        if (i == a.size() - 1 || ccw(p1, a[i], p2)) {
            while(down.size() >= 2 && !ccw(down[down.size()-2], down[down.size()-1], a[i]))
                down.pop_back();
            down.push_back(a[i]);
        }
    }

    a.clear();
    for (int i = 0; i < (int)up.size(); i++)
        a.push_back(up[i]);
    for (int i = down.size() - 2; i > 0; i--)
        a.push_back(down[i]);
}

///////////////////////////////////////////////////

//Get triangle area
double triangle_area(Triangle &t){
	return abs((t.b.x - t.a.x)*(t.c.y - t.a.y) - (t.b.y - t.a.y)*(t.c.x - t.a.x))/2.0;
}

///////////////////////////////////////////////////

//Test if a point is inside a triangle
bool inside_triangle(pt p, Triangle &t){
	Triangle t1, t2, t3;
	t1.a = t.a;
	t1.b = t.b;
	t1.c = p;
	t2.a = t.a;
	t2.b = t.c;
	t2.c = p;
	t3.a = t.b;
	t3.b = t.c;
	t3.c = p;

	return Equal(triangle_area(t1) + triangle_area(t2) + triangle_area(t3), triangle_area(t));
}

///////////////////////////////////////////////////

//Get the convex hull area (useful to test if a point is inside the convex hull)
double convex_hull_area(vector<pt> convex_hull){
	double area = 0.0;

	for(int i=1;i<(int)convex_hull.size()-1;i++){
		Triangle t;
		t.a = convex_hull[0];
		t.b = convex_hull[i];
		t.c = convex_hull[i+1];
		area += triangle_area(t);
	}

	return area; 
}

///////////////////////////////////////////////////

//Test if a point is inside the convex hull
bool inside_convex_hull(pt p, vector<pt> &convex_hull){
	double area_point = 0.0;
	Triangle t;

	for(int i=0;i<(int)convex_hull.size() - 1;i++){	
		t.a = p;
		t.b = convex_hull[i];
		t.c = convex_hull[i+1];
		area_point += triangle_area(t);
	}

	t.a = p;
	t.b = convex_hull[convex_hull.size() - 1];
	t.c = convex_hull[0];
	area_point += triangle_area(t);

	return Equal(area_point, convex_hull_area(convex_hull));
}

///////////////////////////////////////////////////

//Test if point p see edge ab in the convex hull
bool point_see_edge(pt p, pt a, pt b, vector<pt> &convex_hull){
	pt m;
	m.x = (a.x + b.x)/2;
	m.y = (a.y + b.y)/2;
	int ch_size = convex_hull.size();
	int count = 0;

	for(int i=0;i<ch_size;i++){
		pt intersect;
		if(segment_intersect_segment(convex_hull[i], convex_hull[(i+1)%ch_size], p, m, intersect)) count++;
	}

	return count == 1;
}

///////////////////////////////////////////////////

//Get the initial triangle
//If param = 1, search the empty triangle with minimum area
//If param = 2, search the empty triangle with maximum area
vector<pt> get_init_sol(vector<pt> &points, vector<bool> &marked, int param=1){
	if(param == 1){
		double min_area = -1.0;
		int ind_choosen_a, ind_choosen_b, ind_choosen_c;

		for(int i=0;i<(int)points.size();i++){
			for(int j=i+1;j<(int)points.size();j++){
				for(int z=j+1;z<(int)points.size();z++){
					Triangle t;
					t.a = points[i];
					t.b = points[j];
					t.c = points[z];

					bool empty = true;

					for(int k=0;k<(int)points.size();k++){
						if(k != i && k != j && k != z){
							if(inside_triangle(points[k], t)){
								empty = false;
								break;
							}
						}
					}

					if(empty){
						if(min_area == -1.0 || (triangle_area(t) < min_area && !Equal(triangle_area(t), 0.0))){
							min_area = triangle_area(t);
							ind_choosen_a = i;
							ind_choosen_b = j;
							ind_choosen_c = z;
						}
					}
				}
			}
		}

		vector<pt> init_sol;
		init_sol.push_back(points[ind_choosen_a]);
		init_sol.push_back(points[ind_choosen_b]);
		init_sol.push_back(points[ind_choosen_c]);
		marked[ind_choosen_a] = true;
		marked[ind_choosen_b] = true;
		marked[ind_choosen_c] = true;

		return init_sol;
	}
	else if(param == 2){
		double max_area = -1.0;
		int ind_choosen_a, ind_choosen_b, ind_choosen_c;

		for(int i=0;i<(int)points.size();i++){
			for(int j=i+1;j<(int)points.size();j++){
				for(int z=j+1;z<(int)points.size();z++){
					Triangle t;
					t.a = points[i];
					t.b = points[j];
					t.c = points[z];

					bool empty = true;

					for(int k=0;k<(int)points.size();k++){
						if(k != i && k != j && k != z){
							if(inside_triangle(points[k], t)){
								empty = false;
								break;
							}
						}
					}

					if(empty){
						if(triangle_area(t) > max_area){
							max_area = triangle_area(t);
							ind_choosen_a = i;
							ind_choosen_b = j;
							ind_choosen_c = z;
						}
					}
				}
			}
		}

		vector<pt> init_sol;
		init_sol.push_back(points[ind_choosen_a]);
		init_sol.push_back(points[ind_choosen_b]);
		init_sol.push_back(points[ind_choosen_c]);
		marked[ind_choosen_a] = true;
		marked[ind_choosen_b] = true;
		marked[ind_choosen_c] = true;
		
		return init_sol;
	}
}

///////////////////////////////////////////////////

//Renew the current polygonalization by adding one more point to it
//If param = 1, we are searching the polygonalization with minimum area
//If param = 2, we are searching the polygonalization with maximum area
void get_one_step(vector<pt> &cur_sol, vector<pt> &points, vector<bool> &marked, int param=1){
	int cur_sol_size = cur_sol.size();

	if(param == 1){
		double min_area = -1.0;
		int ind_cur_sol, ind_points;

		for(int i=0;i<(int)points.size();i++){
			if(!marked[i]){
				for(int j=0;j<cur_sol_size;j++){
					if(point_see_edge(points[i], cur_sol[j], cur_sol[(j+1)%cur_sol_size], cur_sol)){
						vector<pt> guess_sol;

						for(int z=0;z<=j;z++) guess_sol.push_back(cur_sol[z]);
						guess_sol.push_back(points[i]);
						for(int z=j+1;z<cur_sol_size;z++) guess_sol.push_back(cur_sol[z]);

						vector<pt> ch_guess_sol;
						for(int z=0;z<(int)guess_sol.size();z++) ch_guess_sol.push_back(guess_sol[z]);
						convex_hull(ch_guess_sol);

						bool is_valid = true;

						for(int z=0;z<(int)points.size();z++){
							if(z != i && !marked[z]){
								if(inside_convex_hull(points[z], ch_guess_sol)){
									is_valid=false;
									break;
								}
							}
						}

						if(is_valid){
							Triangle t;
							t.a = cur_sol[j];
							t.b = cur_sol[(j+1)%cur_sol_size];
							t.c = points[i];

							if(min_area == -1.0 || triangle_area(t) < min_area){
								min_area = triangle_area(t);
								ind_cur_sol = j;
								ind_points = i;
							}
						}
					}
				}
			}
		}

		marked[ind_points] = true;
		vector<pt> new_sol;

		for(int i=0;i<=ind_cur_sol;i++) new_sol.push_back(cur_sol[i]);
		new_sol.push_back(points[ind_points]);
		for(int i=ind_cur_sol+1;i<cur_sol_size;i++) new_sol.push_back(cur_sol[i]);

		cur_sol.clear();
		for(int i=0;i<(int)new_sol.size();i++) cur_sol.push_back(new_sol[i]);
	}
	else if(param == 2){
		double max_area = -1.0;
		int ind_cur_sol, ind_points;

		for(int i=0;i<(int)points.size();i++){
			if(!marked[i]){
				for(int j=0;j<cur_sol_size;j++){
					if(point_see_edge(points[i], cur_sol[j], cur_sol[(j+1)%cur_sol_size], cur_sol)){
						vector<pt> guess_sol;

						for(int z=0;z<=j;z++) guess_sol.push_back(cur_sol[z]);
						guess_sol.push_back(points[i]);
						for(int z=j+1;z<cur_sol_size;z++) guess_sol.push_back(cur_sol[z]);

						vector<pt> ch_guess_sol;
						for(int z=0;z<(int)guess_sol.size();z++) ch_guess_sol.push_back(guess_sol[z]);
						convex_hull(ch_guess_sol);

						bool is_valid = true;

						for(int z=0;z<(int)points.size();z++){
							if(z != i && !marked[z]){
								if(inside_convex_hull(points[z], ch_guess_sol)){
									is_valid=false;
									break;
								}
							}
						}

						if(is_valid){
							Triangle t;
							t.a = cur_sol[j];
							t.b = cur_sol[(j+1)%cur_sol_size];
							t.c = points[i];

							if(triangle_area(t) > max_area){
								max_area = triangle_area(t);
								ind_cur_sol = j;
								ind_points = i;
							}
						}
					}
				}
			}
		}

		marked[ind_points] = true;
		vector<pt> new_sol;

		for(int i=0;i<=ind_cur_sol;i++) new_sol.push_back(cur_sol[i]);
		new_sol.push_back(points[ind_points]);
		for(int i=ind_cur_sol+1;i<cur_sol_size;i++) new_sol.push_back(cur_sol[i]);

		cur_sol.clear();
		for(int i=0;i<(int)new_sol.size();i++) cur_sol.push_back(new_sol[i]);
	}
}

///////////////////////////////////////////////////

//Get the solution
//If param = 1, then gets the polygon with minimum area
//If param = 2, then gets the polygon with maximum area
vector<pt> solve(vector<pt> &points, int param=1){
	vector<bool> marked(points.size(), false);
	vector<pt> init_sol = get_init_sol(points, marked, param);

	for(int i=0;i<(int)points.size() - 3;i++) 
		get_one_step(init_sol, points, marked, param);

	return init_sol;
}

///////////////////////////////////////////////////

//Get random double in interval [lower_bound, upper_bound]
double random_double(double lower_bound, double upper_bound){
   	random_device rd;
    mt19937 gen(rd()); 
    uniform_real_distribution<double> dis(lower_bound, upper_bound);
    return dis(gen);
}

///////////////////////////////////////////////////

//Get n points around a circle with radio r
vector<pt> random_circle_points(int n, double r){
	vector<pt> points;

	for(int i=1;i<=n;i++){
		double theta = random_double(0.0, 2*PI);
		pt p;
		p.x = r*cos(theta);
		p.y = r*sin(theta);
		points.push_back(p);
	}

	return points;
}

///////////////////////////////////////////////////

//Get a vector with time needed to solve instances from size n_lower to n_high when jumping step by step
//If param = 1 always search the polygon with minimum area
//If param = 2 always search the polygon with maximum area
pair<vector<double>, vector<double>> get_performance(int n_low=10, int n_high=100, int step=10, int param=1){
	vector<double> x, y;
	double r = 10.0;

	for(int i=n_low;i<=n_high;i+=step){
		x.push_back(1.0*i);
		auto start = std::chrono::high_resolution_clock::now();
		vector<pt> points = random_circle_points(i, r);
		vector<pt> sol = solve(points, param);
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		y.push_back(elapsed.count());
	}

	return {x, y};
}

///////////////////////////////////////////////////

//Plot the time needed to solve instances from size n_lower to n_high when jumping step by step
//If param = 1 always search the polygon with minimum area
//If param = 2 always search the polygon with maximum area
void plot_performance(int n_low=10, int n_high=100, int step=10, int param=1){
	pair<vector<double>, vector<double>> xy = get_performance(n_low, n_high, step, param);
	vector<double> x = xy.first;
	vector<double> y = xy.second;
	plt::figure_size(1200, 780);
	plt::named_plot("Time (in seconds)", x, y);
	plt::title("Performance evaluation");
	string name = "plot";
	if(param == 1) name += "_minimum.png";
	else if(param == 2) name += "_maximum.png";
	plt::save(name);
}

int main(){
	//Uncomment to plot the performance
	//plot_performance(10, 300, 10, 1);
	//plot_performance(10, 300, 10, 2);

	//Example with eight points
	pt p1, p2, p3, p4, p5, p6, p7, p8;
	p1.x = 2.0;
	p1.y = 0.0;
	p2.x = 3.0;
	p2.y = 5.0;
	p3.x = 0.0;
	p3.y = 3.0;
	p4.x = -3.0;
	p4.y = 5.0;
	p5.x = -2.0;
	p5.y = 0.0;
	p6.x = -3.0;
	p6.y = -5.0;
	p7.x = 0.0;
	p7.y = -3.0;
	p8.x = 3.0;
	p8.y = -5.0;

	vector<pt> points;
	points.push_back(p1);
	points.push_back(p2);
	points.push_back(p3);
	points.push_back(p4);
	points.push_back(p5);
	points.push_back(p6);
	points.push_back(p7);
	points.push_back(p8);

	vector<pt> sol_min = solve(points, 1);
	vector<pt> sol_max = solve(points, 2);

	cout << "Solution for the minimum search" << endl;
	for(int i=0;i<(int)sol_min.size();i++) cout << sol_min[i].x << " " << sol_min[i].y << endl;

	cout << "Solution for the maximum search" << endl;
	for(int i=0;i<(int)sol_max.size();i++) cout << sol_max[i].x << " " << sol_max[i].y << endl; 
	return 0;
}