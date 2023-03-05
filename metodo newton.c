// Programa que resuelve sistemas de ecuaciones no lineales usando el método de Neewton//

#include <stdio.h>
#include <math.h>

int main(){
	int opcion,o;
    
	//Mostrar los sitemas y elegir uno//
	do{
		system("cls");
		printf("\t\t\t\tMetodo de Newton para sistemas de ecuaciones no lineales");
		printf("INTEGRANTES DEL EQUIPO:\n Casas Ibañez Emanuel Alexis\n Flores Yanes Francisco Miguel \n Hernández Esquivel Víctor Manuel \n Hernández González Armando \n ");
		printf("\nElige uno de los siguientes sistemas de ecuaciones para resolver\n\n");
		printf("Sistema de Ecuaciones 1\nf1(x, y)= x^2 + xy - 10 = 0\nf2(x, y)= y + 3x^2 - 50 = 0");
		printf("\n\nSistema de Ecuaciones 2\nf1(x, y)= x^2 + y^2 - 9 = 0\nf2(x, y)= -e^x - 2y - 3 = 0");
		printf("\n\nSistema de Ecuaciones 3\nf1(x, y, z)= 2x^2 - 4x + y^2 + 3z^2 + 6z + 2 = 0\nf2(x, y, z)= x^2 + y^2 -2y + 2z^2 -5 = 0\nf3(x, y, z)= 3x^2 - 12x + y^2 -3z^2 + 8= 0");
		printf("\n\nSistema de Ecuaciones 4\nf1(x, y, z)= x^2 + 4x + y^2 = 0\nf2(x, y, z)= x^2 -x -12y + 1 = 0\nf3(x, y, z)= 3x^2 -12x + y^2 -3z^2 + 8 = 0");
		printf("\n\n 5: Salir");
		printf("\n\nElige una opcion (1,...,5): ");
		scanf("%i",&opcion);
		if(opcion == 5) printf("\nAdios");
		else if(opcion > 5) printf("No elegiste una opcion correcta");
		else{
			newton(opcion);
			printf("\n\nDeseas probar con otros puntos iniciales?\nSi:1 No:2\n.");
			scanf("%i",&o);
			if(o == 1){
				newton(opcion);
			}
			printf("\n\nDeseas resolver otro sitema?\nSi:1 No:2\n.");
			scanf("%i",&o);
			if(o == 2) opcion = 5;
		}
	}while(opcion != 5);
}


void newton(int op){
	int ite,i;
	double x,y,z,det,x1,y1,z1,x2,y2,z2,xt,yt,mayor,mayor1;
	double e = 2.71828183;
	double tol;
	double jac1[5][5];
	double jac[5][5];
	double fun[4][3];
	
	if(op <= 2){
		printf("Determinar puntos iniciales (x,y)");
		printf("\nx: ");
		scanf("%lf",&x);
		printf("y: ");
		scanf("%lf",&y);
		printf("Numero de iteraciones que se desean realizar: ");
		scanf("%i",&ite);
		printf("Tolerancia: ");
		scanf("%lf",&tol);
		for(i=1; i<=ite; i++){
			if(op == 1){
				jac[0][0] = (2*x) + y;
				jac[0][1] = x;
				jac[1][0] = 3 * y*y;
				jac[1][1] = 6*x*y + 1;
				fun[0][0] = (x*x) + x*y - 10;
				fun[1][0] = y + 3*x*y*y - 50;
		    }else if(op == 2){
	    		jac[0][0] = 2*x;
	    		jac[0][1] = 2*y;
	    		jac[1][0] = pow(-e,x);
	    		jac[1][1] = -2;
	    		fun[0][0] = pow(x,2) + pow(y,2) -9;
	    		fun[1][0] = pow(-e,x) - 2*y - 3;
			}
			printf("\n\nJ(X) = \t\t\t\t\t F(X) = ");
			printf("\n|%f\t%f|\t\t|%f|\n|%f\t%f|\t\t|%f|",jac[0][0],jac[0][1],fun[0][0],jac[1][0],jac[1][1],fun[1][0]);
			
			
			det = ((jac[0][0]*jac[1][1]-jac[1][0]*jac[0][1]));
			jac[2][2] = jac[0][0];
			jac[0][0] = (1/det)*jac[1][1];
			jac[0][1] = (1/det)*(-(jac[0][1]));
			jac[1][0] = (1/det)*(-(jac[1][0]));
			jac[1][1] = (1/det)*jac[2][2];
			
		
			x1 =  (jac[0][0]*fun[0][0]) + (jac[0][1]*fun[1][0]);
			y1 =  (jac[1][0]*fun[0][0]) + (jac[1][1]*fun[1][0]);
			x2 = x;
			y2 = y;
			x = x - x1;
			y = y - y1;
			
			
			printf("\n\nIteracion %i \tx = %f\t\ty = %f",i,x,y);
			
			xt = fabs(x2) - fabs(x);
			yt = fabs(y2) - fabs(y);
			if(xt > yt){
				mayor = xt;
			}else{
				mayor = yt;
			}
			if(fabs(x) > fabs(y)){
				mayor1 = fabs(x);
			}else{
				mayor1 = fabs(y);
			}
			
			
			
			tol = mayor/mayor1;
			printf("\n\n El error es %f",tol);
			
			 
			
		}
	}else{
		printf("Determinar puntos iniciales (x,y,z)");
		printf("\nx: ");
		scanf("%lf",&x);
		printf("y: ");
		scanf("%lf",&y);
		printf("z: ");
		scanf("%lf",&z);
		printf("Numero de iteraciones que se desean realizar: ");
		scanf("%i",&ite);
		printf("Tolerancia: ");
		scanf("%lf",&tol);
		for(i = 1; i <= ite; i++){
		
			if(op==3){
				jac[0][0] = (4*x) - 4;
				jac[0][1] = 2*y;
				jac[0][2] = (6*z) + 6;
				jac[1][0] = 2*x;
				jac[1][1] = (2*y) - 2;
				jac[1][2] = 4*z;
				jac[2][0] = (6*x) - 12;
				jac[2][1] = 2*y;
				jac[2][2] = -6*z;
				fun[0][0] = (2*pow(x,2)) - 4*x + pow(y,2) + (3*pow(z,2) + 6*z + 2);
				fun[1][0] = pow(x,2) + pow(y,2) - (2*y) + 2*pow(z,2) - 5;
				fun[2][0] = (3*pow(x,2)) - 12*x + pow(y,2) - (3*pow(z,2)) + 8;
			}else if(op == 4){
				jac[0][0] = (2*x) - 4;
				jac[0][1] = 2*y;
				jac[0][2] = 0;
				jac[1][0] = (2*x) - 1;
				jac[1][1] = -12;
				jac[1][2] = 0;
				jac[2][0] = (6*x) - 12;
				jac[2][1] = 2*y;
				jac[2][2] = -6*z;
				fun[0][0] = (pow(x,2)) - 4*x + (pow(y,2));
				fun[1][0] = pow(x,2) - x - 12*y + 1;
				fun[2][0] = (3*pow(x,2)) - 12*x + pow(y,2) - (3*pow(z,2)) + 8;
			}
			printf("\n\nJ(X) = \t\t\t\t\t\t\t F(X) = ");
			printf("\n|%f\t%f\t%f|\t\t|%f|\n|%f\t%f\t%f|\t\t|%f|\n|%f\t%f\t%f|\t\t|%f|",jac[0][0],jac[0][1],jac[0][2],fun[0][0],jac[1][0],jac[1][1],jac[1][2],fun[1][0],jac[2][0],jac[2][1],jac[2][2],fun[2][0]);
			
			det = (jac[0][0]*(jac[1][1]*jac[2][2]-jac[2][1]*jac[1][2])) - (jac[0][1]*(jac[1][0]*jac[2][2] - jac[2][0]*jac[1][2])) + (jac[0][2]*(jac[1][0]*jac[2][1] - jac[2][0]*jac[1][1]));
			jac1[0][0] = (jac[1][1]*jac[2][2]) - (jac[2][1]*jac[1][2]);
			jac1[0][1] = (jac[0][1]*jac[2][2]) - (jac[2][1]*jac[0][2]);
			jac1[0][2] = (jac[0][1]*jac[1][2]) - (jac[1][1]*jac[0][2]);
			jac1[1][0] = (jac[1][0]*jac[2][2]) - (jac[2][0]*jac[1][2]);
			jac1[1][1] = (jac[0][0]*jac[2][2]) - (jac[2][0]*jac[0][2]);
			jac1[1][2] = (jac[0][0]*jac[1][2]) - (jac[1][0]*jac[0][2]);
			jac1[2][0] = (jac[1][0]*jac[2][1]) - (jac[2][0]*jac[1][1]);
			jac1[2][1] = (jac[0][0]*jac[2][1]) - (jac[2][0]*jac[0][1]);
			jac1[2][2] = (jac[0][0]*jac[1][1]) - (jac[1][0]*jac[0][1]);
			
			jac[0][0] = (1/det)  * jac1[0][0];
			jac[0][1] = -(1/det) * jac1[0][1];
			jac[0][2] = (1/det)  * jac1[0][2];
			jac[1][0] = -(1/det) * jac1[1][0];
			jac[1][1] = (1/det)  * jac1[1][1];
			jac[1][2] = -(1/det) * jac1[1][2];
			jac[2][0] = (1/det)  * jac1[2][0];
			jac[2][1] = -(1/det) * jac1[2][1];
			jac[2][2] = (1/det)  * jac1[2][2];
			
			
			x1 = (jac[0][0]*fun[0][0]) + (jac[0][1]*fun[1][0]) + (jac[0][2]*fun[2][0]);
			y1 = (jac[1][0]*fun[0][0]) + (jac[1][1]*fun[1][0]) + (jac[1][2]*fun[2][0]);
			z1 = (jac[2][0]*fun[0][0]) + (jac[2][1]*fun[1][0]) + (jac[2][2]*fun[2][0]);
			
			x = x - x1;
			y = y - y1;
			z = z - z1;
			printf("\n\nIteracion %i \tx = %f\t\ty = %f\t\tz = %f",i,x,y,z);
			
		}
		
		
	}
}

	
	














