/*
TODO:
Change vector pairs to a class that returns magnitude
Convert Everything to AU
Figure out way to graph 
*/

#include <iostream>
#include <cmath>
#include <SFML/Graphics.hpp>
using namespace std;

const int WinX = 1920;
const int WinY = 1080;

//Gravitaional Constant in m^3 kg^-1 s^-2
const double G = 6.6730*pow(10,20);
//Mass of Sun in kg
const double Ms = 1.989*pow(10,30);
//Initial Velocity of Sattelite in m/s
const double Vi = 29722;
//Initial Position of Sattelite in m (x, 0)
const double Xi = 150000000;
//Step size in seconds
const int h = 86400;
//# of steps
const int t = 1000;

//radial acceleration function (Newtonian)
double a(double r){
	return G*Ms/(r*r);
}
pair<double,double> vectorizeA(double a,pair<double,double> v){
	//perependicular v
	pair<double,double> vp;
	vp.first = v.second*a / pow((pow(v.first,2)+pow(v.second,2)),.5);
	vp.second = -1*v.first*a / pow((pow(v.first,2)+pow(v.second,2)),.5);
	return vp;
}
pair<double,double> vectorAdd(pair<double,double> i, pair<double,double> j){
	pair<double,double> k (i.first + j.first,i.second + j.second);
	return k;

}
pair<double,double> vectorAdd4(pair<double,double> i, pair<double,double> j, pair<double,double> k, pair<double,double> l){
	pair<double,double> m (i.first + j.first +k.first + l.first  ,i.second + j.second +k.second + l.second );
	return m;
}
pair<double,double> vectorMultiply(pair<double,double>i, double j){
	pair<double,double> k (i.first*j,i.second*j);
	return k;
}
//rectangular formalization
pair <double,double> v (0, Vi);
pair <double,double> pos (Xi,0);

bool runge = false;
//sun at 0,0, assuming Msun >> Msat
int main(){
	//SFML Window
	 // create the window
    sf::RenderWindow window(sf::VideoMode(WinX, WinY), "My window");

    // run the program as long as the window is open
	int time = 0;
	int rad = 100;
	sf::CircleShape Sun(rad);
	Sun.setFillColor(sf::Color(252, 212, 64));
	Sun.setPosition(sf::Vector2f(WinX/2,WinY/2));
	window.clear(sf::Color::Black);
    while (window.isOpen())
    {
		\ss
        // check all the window's events that were triggered since the last iteration of the loop
        sf::Event event;
        while (window.pollEvent(event))
        {
            // "close requested" event: we close the window
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // clear the window with black color
        window.clear(sf::Color::Black);
		
		
        // draw everything here...
        window.draw(Sun);
        // end the current frame
        window.display();
    }

	
	
	
	
	
	//Text Output
	freopen("output.txt", "w", stdout);
	//euler implementation
	if(!runge){
	for (int z = 0;z<t;z++){
		double r1 = pow((pow(pos.first,2)+pow(pos.second,2)),.5);
		//define r2 when doing rk4 or something
		v = vectorAdd(v,vectorMultiply(vectorizeA(a(r1),v),h));
		pos = vectorAdd(pos,vectorMultiply(v,h));
		cout << pos.first << " "<< pos.second <<endl;
	}}/*
	}
	else{
		double r1 = pow((pow(pos.first,2)+pow(pos.second,2)),.5);
		//define r2 when doing rk4 or something
		v = vectorAdd(v,vectorMultiply(vectorizeA(a(r1),v),h));
		pos = vectorAdd(pos,vectorMultiply(v,h));
		
	}*/
	return 0;
}