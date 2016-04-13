import org.jfree.ui.RefineryUtilities; 
import java.util.Arrays;


public class TheoryModel{

	////  ** Arbitrary Values Taken Here ** ////

	// FOUNDATIONAL THEORY MODEL 
	public static final int T = 20 ;						// number of years in a generation

	static float K_0 = 136.0f;								// starting physical capital 
	static float KN_0 = 50.0f ;							// starting knowledge capital 
	static float Y_0 = 47.36f;								// starting GDP

	static float alpha = 0.35f;
	static float beta = 0.55f;								//  0 < alpha, beta < 1

	static float s = 0.3f;									// constant fraction of output that is saved and invested
	static float delta = 0.0005f;							// exogenously fixed capital depreciation rate
	static float n = 0.0005f;								// exogenous growth rate of cohort
	static float LF_0 = 12.0f;								//assumed Labour Force at start of base data 
	public static float B_0 = 3.0f;						//assumed Cohort size at the start of the base data
	static int theta_0 = 3;

	// THE ROLE OF THE GOVERNMENT
	static float eta = 0.2f;								// constant fraction of income y\(t\) that is invested in the RnD
	static float phi = 5.0f;								// phi > 0, conversion rate of education to cohort specific absorpition capacity 

	static float u_upper_0 = 0.06f;
	static float u_lower_0 = 0.01f;							// u > 0; u_upper > u_lower; u_lower != 0;

	static float psi_const = 0.33f;
	static float gamma_const = 0.33f;
	static int tau = 250;									// secondary time variable
	static float rho = 0.0001f; 							// time discount rate
	static float rd_0 = 0.0125f;							// investment in rnd for startinf year

	////  ** Arbitrary Values Taken Here ** ////

// FOUNDATIONAL THEORY MODEL 

	//  Cobb Douglas type Production Function variables
	public static float[] Y_t = new float[1000];			// GDP of Economy at time t
	public static float[] K_t = new float[1000];			// Physical Capital at time t
	public static float[] KN_t = new float[1000];			// Knowledge Capital at time t
	public static float[] G1_t = new float[1000];			// Human Capital of the first generation at time t
	public static float[] G2_t = new float[1000];			// Human Capital of the second generation at time t
	public static float[] G3_t = new float[1000];			// Human Capital of the third generation at time t


// EVOLUTIONARY VARIABLES

	public static float[] B_t = new float[1000];			// new Cohort born at each time t ( number of people born at time t )
	public static float[] L_t = new float[1000];			// size of the population at time 	

	public static float[] y_t = new float[1000];			// Per Capita GDP of Economy at time t
	public static float[] k_t = new float[1000];			// Physical Capital per worker at time t
	public static float[] kn_t = new float[1000];			// Per Capita Knowledge Capital at time t

	public static int[] theta_t = new int[1000];

	public static float[] LF_t = new float[1000];


// THE ROLE OF THE GOVERNMENT
	public static float[] c_t = new float[1000];			//cost of implementation of unit technology
	public static float[] i_t = new float[1000];			// intermediate good
	public static float[] rd_t = new float[1000];			// resources devoted to RnD sector

	public static float[] h_t = new float[1000];			// quality adjusted human capital of an individual (human capital of a representative worker)
	public static float[] a_t = new float[1000];			// absorption capacity of the cohort
	public static float[] H_t = new float[1000];			// cohort specific human capital 



	public static float calculate_v_t(int year, float[] L_t, float[] G1_t, float psi_const, float n, float delta,
								 float s, float[] G2_t, float gamma_const, float[] G3_t, float K_t, 
								 float alpha, float beta, int theta_iter, float rho){
		float KN;
		float K = K_t;
		float Y=0.0f;
		float sum=0.0f;
		int k;
		for (k=0;k<=theta_iter;k++){
			KN = (float)(Math.pow(G1_t[year+k],psi_const)*Math.pow(G2_t[year+k],gamma_const)*Math.pow(G3_t[year+k],1-psi_const-gamma_const));
			Y = (float)(Math.pow(K,alpha)*Math.pow(KN,beta));
			sum = sum + (float)(Math.exp(k*rho)*Y/L_t[year+k]);
			K = K + s*Y - (n+delta)*K;
		}

		float temp = Y/L_t[year+k-1];
		for (k=theta_iter+1; k<100;k++){
			sum = sum + (float)(Math.exp(k*rho)*temp);
		}

		return sum;
	}



	public static void main(String[] args){

		// FOUNDATIONAL THEORY MODEL 
		int i,j,k;
		int edu = 0;
		float add;
		K_t[0] = K_0;					// starting physical capital
		KN_t[0] = KN_0;					// starting knowledge capital
		Y_t[0] = Y_0;					// starting GDP
		B_t[0] = B_0;					// starting cohort ( also starting population )
		L_t[0] = B_0+LF_0;					// starting population
		k_t[0] = K_0/B_0;				// per capita starting physical capital
		kn_t[0] = KN_0/B_0;				// per capita starting knowledge capital 
		y_t[0] = Y_0/B_0;				// per capita GDP
		rd_t[0] = rd_0;		 			// starting investment in rnd
		c_t[0] = y_t[0]*10;				//starting cost of technology implementation

		// theta(t)	:: 0 - 19 yrs
		for (i=0;i<T;i++){
			theta_t[i] = theta_0;
		}

		// setting assumed Labour Force at start of base data until people 
		// start joining after starting theta_t(0) yrs
		// setting assumed starting Knowledge Capital
		for(i=0;i<3*T;i++){
			if (i<T){
				LF_t[i] = LF_t[i] + LF_0;
				G1_t[i] = G1_t[i] + (float)(Math.pow(KN_0,0.33));
				G2_t[i] = G2_t[i] + (float)(Math.pow(KN_0,0.33));
				G3_t[i] = G3_t[i] + (float)(Math.pow(KN_0,0.33));
			}
			else if (i>=T && i<2*T){
				LF_t[i] = LF_t[i] + LF_0*(2/3);
				G2_t[i] = G2_t[i] + (float)(Math.pow(KN_0,0.33));
				G3_t[i] = G3_t[i] + (float)(Math.pow(KN_0,0.33));
			}
			else{
				LF_t[i] = LF_t[i] + LF_0/3;
				G3_t[i] = G3_t[i] + (float)(Math.pow(KN_0,0.33));
			}
		}

		// KN(t), K(t), Y(t), L(t), B(t), kn(t), k(t), y(t) ::	0 - 80 yrs
		for (i = 0; i < 4*T; i++){
			// theta(t)	:: 20 - 79 yrs
			if (i>=T){
				if (i%T == 0){
					// education value is set once every 20 years(generation lenght) and not every year for the period of starting data
					edu = (int)(Math.random()*4)  + 3;			
				}
				theta_t[i] = edu;
			}


			// THE ROLE OF THE GOVERNMENT
			// c(t), rd(t), i(t), a(t), h(t), H(t) :: 0 - 80 yrs
			if (i>0){
				c_t[i] = c_t[i-1] + ((float)Math.random()*0.015f)*c_t[i-1];
				rd_t[i] = eta * y_t[i-1];
			}
			i_t[i] = u_lower_0*rd_t[i];									//  (7)
			a_t[i] = (float)Math.exp(phi*theta_t[i]);
			h_t[i] = a_t[i] * i_t[i];									//   (9)
			H_t[i] = h_t[i] * B_t[i];	


			for(k=theta_t[i];k<3*T;k++){
				LF_t[i+k] = LF_t[i+k] + B_t[i];
				if (k<T){
					G1_t[i+k] = G1_t[i+k] + H_t[i];
				}
				else if (k>=T && k<2*T){
					G2_t[i+k] = G2_t[i+k] + H_t[i];
				}
				else{
					G3_t[i+k] = G3_t[i+k] + H_t[i];
				}
			}
			
			
			//  Cobb Douglas type Production Function
			K_t[i+1] = K_t[i] + s*Y_t[i] - (n+delta)*K_t[i];
			KN_t[i] = (float)(Math.pow(G1_t[i],psi_const)*Math.pow(G2_t[i],gamma_const)*Math.pow(G3_t[i],1-psi_const-gamma_const));
			Y_t[i] = (float)Math.pow( K_t[i] , alpha ) * (float)Math.pow( KN_t[i], beta);		//  (1)

			//population and cohort growth
			B_t[i+1] = (B_t[i]*(float)Math.exp(n*1));
			L_t[i+1] = L_t[i] + B_t[i+1];

			k_t[i+1] = K_t[i+1]/L_t[i+1];
			kn_t[i] = KN_t[i]/L_t[i];
			//  Per Capita Cobb Douglas Production Function
			y_t[i] = (float)Math.pow( k_t[i] , alpha ) *(float)Math.pow( kn_t[i], beta);			//  (2)	
			

		}


		// int index=0;
		// for(float s : G1_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		}
		// index=0;
		// for(float s : G2_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }
		// }
		// index=0;
		// for(float s : G3_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		}
		// index=0;
		// for(float s : y_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		
		// }

		// index=0;
		// for(float s : L_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		
		// }
		// index=0;
		// for(float s : i_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		
		// }
		// index=0;
		// for(float s : a_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		
		// }
		// index=0;
		// for(float s : h_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		
		// }
		// index=0;
		// for(float s : theta_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		
		// }
		// index=0;
		// for(float s : k_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		
		// }
		// index=0;
		// for(float s : kn_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		
		// }
		// index=0;
		// for(float s : y_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		
		// }

// ALL Evolutionary BASE DATA STORED ARE APPROPRIATELY GENERATED BY NOW FROM T = 0 TO T = 80 YRS
	// primary time variable
		int t = 4*T;				

	// since B(t) and L(t) follow straight forward trend I am calculating it here separately for period 80 - 680 yrs
		for (j=0;j<20*T;j++){
			B_t[t+j] = B_t[t+j-1] * (float)Math.exp(n*1);
			L_t[t+j] = L_t[t+j-1] + B_t[t+j] - B_t[j];
		}


//THE MAIN MODEL STARTS WORKING HERE ONWARDS

		boolean techJump =  false;		// indicator of whether technology jump has happened or not
		float w_t=0.0f, v_t = 0.0f, sum4=0.0f, max_w_t=0.0f;
		int theta_iter=0, max_theta=0;

		for(j=0;j<3*T;j++){

			max_w_t = -100.0f;
			max_theta = 0;

			c_t[t+j] = c_t[t+j-1] + ((float)Math.random()*0.015f)*c_t[t+j-1];
			rd_t[t+j] = eta * y_t[t+j-1];

			if (!techJump){		
			//jump has not happened yet
				for (theta_iter = 0; theta_iter<T; theta_iter++){
				// for calculating with u_lower
					i_t[t+j] = u_lower_0*rd_t[t+j];									//  (7)
					a_t[t+j] = (float)Math.exp(phi*theta_iter);
					h_t[t+j] = a_t[t+j] * i_t[t+j];									//  (9)
					H_t[t+j] = h_t[t+j] * B_t[t+j];

					G1_t[t+j+theta_iter] = G1_t[t+j+theta_iter] + H_t[t+j];
					v_t = calculate_v_t(t+j, L_t, G1_t,psi_const, n, delta, s, G2_t, gamma_const, G3_t, K_t[t+j], 
								 			alpha, beta, theta_iter, rho);
				
					G1_t[t+j+theta_iter] = G1_t[t+j+theta_iter] - H_t[t+j];
					w_t = v_t-0.0001f*theta_iter;
					System.out.println(t+j);
					System.out.println(v_t + " u_lower " + theta_iter );

					if (max_w_t < w_t){
						// System.out.println(w_t + " Yippee lower " +theta_iter);
						max_w_t = w_t;
						max_theta = theta_iter;
					}
				}		
				

				for (theta_iter = 0; theta_iter<T; theta_iter++){
				// for calculating with u_upper
					i_t[t+j] = u_upper_0*rd_t[t+j];									//  (7)
					a_t[t+j] = (float)Math.exp(phi*theta_iter);
					h_t[t+j] = a_t[t+j] * i_t[t+j];									//  (9)
					H_t[t+j] = h_t[t+j] * B_t[t+j];
			
					G1_t[t+j+theta_iter] = G1_t[t+j+theta_iter] + H_t[t+j];
					v_t = calculate_v_t(t+j, L_t, G1_t,psi_const, n, delta, s, G2_t, gamma_const, G3_t, K_t[t+j], 
								 			alpha, beta, theta_iter, rho);
				
					G1_t[t+j+theta_iter] = G1_t[t+j+theta_iter] - H_t[t+j];
					
					w_t = v_t - 0.0001f*theta_iter - c_t[t+j]*(u_upper_0 - u_lower_0);

					System.out.println(v_t + " u_upper " + theta_iter );
					// System.out.println(w_t+"#");
					// System.out.println(theta_iter);

					if (max_w_t < w_t){

						// System.out.println(w_t + " Yippee upper " +theta_iter);
						max_w_t = w_t;
						max_theta = theta_iter;
						techJump = true;
					}
						
				}
			}
			else{
			// the jump has already been made
				for (theta_iter = 0; theta_iter<T; theta_iter++){
					// for calculating with u_upper
					i_t[t+j] = u_upper_0*rd_t[t+j];									//  (7)
					a_t[t+j] = (float)Math.exp(phi * theta_iter);
					h_t[t+j] = a_t[t+j] * i_t[t+j];									//  (9)
					H_t[t+j] = h_t[t+j] * B_t[t+j];
			
					G1_t[t+j+theta_iter] = G1_t[t+j+theta_iter] + H_t[t+j];
					v_t = calculate_v_t(t+j, L_t, G1_t,psi_const, n, delta, s, G2_t, gamma_const, G3_t, K_t[t+j], 
								 			alpha, beta, theta_iter, rho);
				
					G1_t[t+j+theta_iter] = G1_t[t+j+theta_iter] - H_t[t+j];
					w_t = v_t - 0.0001f*theta_iter;

					if (max_w_t < w_t){
						max_w_t = w_t;
						max_theta = theta_iter;
					}
					
				}
			}

			// calculating the parameters based on theta decided
			theta_t[t+j] = max_theta;
			if (!techJump){
				i_t[t+j] = u_lower_0*rd_t[t+j];									//  (7)
			}
			else{
				i_t[t+j] = u_upper_0*rd_t[t+j];	
			}
			a_t[t+j] = (float)Math.exp(phi * theta_t[t+j]);
			h_t[t+j] = a_t[t+j] * i_t[t+j];									//  (9)
			H_t[t+j] = h_t[t+j] * B_t[t+j];

			for(k=theta_t[t+j];k<3*T;k++){
				LF_t[t+j+k] = LF_t[t+j+k] + B_t[t+j];
				if (k<T){
					G1_t[t+j+k] = G1_t[t+j+k] + H_t[t+j];
				}
				else if (k>=T && k<2*T){
					G2_t[t+j+k] = G2_t[t+j+k] + H_t[t+j];
				}
				else{
					G3_t[t+j+k] = G3_t[t+j+k] + H_t[t+j];
				}
			}

			//  Cobb Douglas type Production Function
			K_t[t+j+1] = K_t[t+j] + s*Y_t[t+j] - (n+delta)*K_t[t+j];
			KN_t[t+j] = (float)(Math.pow(G1_t[t+j],psi_const)*Math.pow(G2_t[t+j],gamma_const)*Math.pow(G3_t[t+j],1-psi_const-gamma_const));
			Y_t[t+j] = (float)Math.pow( K_t[t+j] , alpha ) * (float)Math.pow( KN_t[t+j], beta);

			k_t[t+j+1] = K_t[t+j+1]/L_t[t+j+1];
			kn_t[t+j] = KN_t[t+j]/L_t[t+j];
			//  Per Capita Cobb Douglas Production Function
			y_t[t+j] = (float)Math.pow( k_t[t+j] , alpha ) *(float)Math.pow( kn_t[t+j], beta);
			
		}


		// int index=0;
		// for(float s : G1_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		}
		// index=0;
		// for(float s : G2_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }
		// }
		// index=0;
		// for(float s : G3_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		}
		// index=0;
		// for(float s : y_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		
		// }

		int index=0;
		for(float s : L_t){
			if (index<200){
		    	System.out.println(String.valueOf(index++)+": "+s);
		    }		
		}
		// index=0;
		// for(float s : i_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		
		// }
		// index=0;
		// for(float s : a_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		
		// }
		// index=0;
		// for(float s : h_t){
		// 	if (index<200){
		//     	System.out.println(String.valueOf(index++)+": "+s);
		//     }		
		// }
		index=0;
		for(float s : theta_t){
			if (index<200){
		    	System.out.println(String.valueOf(index++)+": "+s);
		    }		
		}
		index=0;
		for(float s : k_t){
			if (index<200){
		    	System.out.println(String.valueOf(index++)+": "+s);
		    }		
		}
		index=0;
		for(float s : kn_t){
			if (index<200){
		    	System.out.println(String.valueOf(index++)+": "+s);
		    }		
		}
		index=0;
		for(float s : y_t){
			if (index<200){
		    	System.out.println(String.valueOf(index++)+": "+s);
		    }		
		}
	
		

		// XYLineChart_AWT_float chart1 = new XYLineChart_AWT_float(" Figure", "Y(t) = K(t)^alpha * KN(t)^beta", "Y(t)", Y_t, 20);
		// chart1.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart1 );          
		// chart1.setVisible( true );

		// XYLineChart_AWT_float chart2 = new XYLineChart_AWT_float(" Figure", "K(t)", "K(t)", K_t, 20);
		// chart2.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart2 );          
		// chart2.setVisible( true );

		// XYLineChart_AWT_float chart3 = new XYLineChart_AWT_float(" Figure", "KN(t)", "KN(t)", KN_t, 20);
		// chart3.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart3 );          
		// chart3.setVisible( true );

		// XYLineChart_AWT_float chart4 = new XYLineChart_AWT_float(" Figure", "B(t)", "B(t)", B_t, 20);
		// chart4.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart4 );          
		// chart4.setVisible( true );

		// XYLineChart_AWT_float chart5 = new XYLineChart_AWT_float(" Figure", "L(t)", "L(t)", L_t, 20);
		// chart5.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart5 );          
		// chart5.setVisible( true );

		// XYLineChart_AWT_float chart6 = new XYLineChart_AWT_float(" Figure", "LF(t)", "LF(t)", LF_t, 20);
		// chart6.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart6 );          
		// chart6.setVisible( true );

		// XYLineChart_AWT_float chart7 = new XYLineChart_AWT_float(" Figure", "k(t)", "k(t)", k_t, 20);
		// chart7.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart7 );          
		// chart7.setVisible( true );

		// XYLineChart_AWT_float chart8 = new XYLineChart_AWT_float(" Figure", "y(t)", "y(t)", y_t, 20);
		// chart8.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart8 );          
		// chart8.setVisible( true );

		// XYLineChart_AWT_float chart9 = new XYLineChart_AWT_float(" Figure", "kn(t)", "kn(t)", kn_t, 20);
		// chart9.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart9 );          
		// chart9.setVisible( true );

		// XYLineChart_AWT_int chart10 = new XYLineChart_AWT_int(" Figure", "theta_t(t)", "theta_t(t)", theta_t, 20);
		// chart10.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart10 );          
		// chart10.setVisible( true );


		// XYLineChart_AWT_float chart16 = new XYLineChart_AWT_float(" Figure", "i(t)", "i(t)", i_t, 20);
		// chart16.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart16 );          
		// chart16.setVisible( true );

		// XYLineChart_AWT_float chart17 = new XYLineChart_AWT_float(" Figure", "a(t)", "a(t)", a_t, 20);
		// chart17.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart17 );          
		// chart17.setVisible( true );

		// XYLineChart_AWT_float chart18 = new XYLineChart_AWT_float(" Figure", "h(t)", "h(t)", h_t, 20);
		// chart18.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart18 );          
		// chart18.setVisible( true );

		// XYLineChart_AWT_float chart19 = new XYLineChart_AWT_float(" Figure", "rd(t)", "rd(t)", rd_t, 20);
		// chart19.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart19 );          
		// chart19.setVisible( true );

		// XYLineChart_AWT_float chart20 = new XYLineChart_AWT_float(" Figure", "H(t)", "H(t)", H_t, 20);
		// chart20.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart20 );          
		// chart20.setVisible( true );

		// XYLineChart_AWT_float chart21 = new XYLineChart_AWT_float(" Figure", "v(t)", "v(t)", v_t, 20);
		// chart21.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart21 );          
		// chart21.setVisible( true );
	}
	
}