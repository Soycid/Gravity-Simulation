/*
TODO:
Change vector pairs to a class that returns magnitude
Convert Everything to AU
Figure out way to graph 
 >add zoom feature
*/

#include <iostream>
#include <cmath>
#include <SFML/Graphics.hpp>
using namespace std;
//one pixel = 300k km
//one time unit = 1 s
const int WinX = 2*1920;
const int WinY = 2*1080;

const double SunX = WinX/2;
const double SunY = WinY/2;
//Gravitaional Constant in pix^3 kg^-1 s^-2
const double G = 2.497413*pow(10,-36);
//Mass of Sun in kg
const double Ms = 1.989*pow(10,30);
//Initial Velocity of Sattelite in pix/s
const double Vi = pow(10,-4);
//Initial Position of Sattelite in pix (x, 0)
const double Xi = 500;
//Step size in seconds
const double h = 86400;
//framerate
const int framerate = 30;
//defines parameters of planet
class Planet {
	public:
		double posX,posY,velX,velY,r,a,ax,ay;
		void euler_step();
};
//defines Euler Method to update position
void Planet::euler_step(){
	r = pow(pow((posX-SunX),2) + pow((posY-SunY),2),.5);
	a = (G*Ms)/(r*r);
	ax = a * (SunX - posX)/r;
	ay = a * (SunY - posY)/r;
	velX += (ax*h);
	velY += (ay*h);
	posX += (velX*h);
	posY += (velY*h);
	
}


bool runge = false;
//assuming Msun >> Msat
int main(){
	freopen("output.txt", "w", stdout);
	//define Earth
	Planet E;
	E.posX = WinX/2 + Xi;
	E.posY = WinY/2;
	E.velX = 0;
	E.velY = 1*Vi;
	
	//SFML Window
	 // create the window
    sf::RenderWindow window(sf::VideoMode(WinX, WinY), "EEgrav");
	window.setFramerateLimit(framerate);
    // run the program as long as the window is open
	int time = 0;
	int rad = 50;
	sf::CircleShape Sun(rad);
	Sun.setFillColor(sf::Color(252, 212, 64));
	Sun.setPosition(sf::Vector2f(SunX,SunY));
	
	sf::CircleShape Earth(rad/10);
	Earth.setFillColor(sf::Color(0, 119, 190));
	Earth.setPosition(sf::Vector2f(SunX + Xi,SunY));
	
	
	
	window.clear(sf::Color::Black);
    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }
		cout << E.posX << " "<< E.posY <<endl;

		E.euler_step();
		int pixX = (int) E.posX;
		int pixY = (int) E.posY;
		Earth.setPosition(sf::Vector2f(pixX,pixY));
        window.clear(sf::Color::Black);		
        window.draw(Sun);
		window.draw(Earth);		
        window.display();
    }

	return 0;
}