import org.jfree.ui.RefineryUtilities; 
import java.util.Arrays;


public class TheoryModel{

	////  ** Arbitrary Value Taken Here ** ////

	// FOUNDATIONAL THEORY MODEL 
	public static final int T = 20 ;										// number of years in a generation

	static float K_0 = 1000.0f;								// starting physical capital 
	static float sum_var2 = 600.0f ;						// starting knowledge capital 
	static float Y_0 = 15.0f;								// starting GDP

	static float alpha = 0.35f;
	static float beta = 0.55f;								//  0 < alpha, beta < 1

	static float s = 0.3f;									// constant fraction of output that is saved and invested
	static float delta = 0.005f;							// exogenously fixed capital depreciation rate
	static float n = 0.005f;									// exogenous growth rate of cohort
	static long LF_0 = 10;									//assumed Labour Force at start of base data 
	public static long B_0 = 10;							//assumed Cohort size at the start of the base data

	// THE ROLE OF THE GOVERNMENT
	static float eta = 0.01f;								// constant fraction of income y\(t\) that is invested in the RnD
	static float phi = 0.068f;								// phi > 0, conversion rate of education to cohort specific absorpition capacity 

	static float udot_0 = 0.05f;								// size of the radical change in technology
	static float u_upper_0 = 0.06f;
	static float u_lower_0 = 0.01f;
	static float u_0 = 0.07f;	

	static float psi_const = 0.33f;
	static float gamma_const = 0.33f;

	static int tau = 250;									// secondary time variable
	static float rho = 0.0001f; 								// discount time
	////  ** Arbitrary Value Taken Here ** ////

// FOUNDATIONAL THEORY MODEL 

	//  Cobb Douglas type Production Function variables
	public static float[] Y_t = new float[1000];				// GDP of Economy at time t
	public static float[] K_t = new float[1000];				// Physical Capital at time t
	public static float[] KN_t = new float[1000];			// Knowledge Capital at time t

// EVOLUTIONARY VARIABLES

	public static long[] B_t = new long[1000];				// new Cohort born at each time t ( number of people born at time t )
	public static long[] L_t = new long[1000];				// size of the population at time 	

	public static float[] y_t = new float[1000];				// Per Capita GDP of Economy at time t
	public static float[] k_t = new float[1000];				// Physical Capital per worker at time t
	public static float[] kn_t = new float[1000];			// Per Capita Knowledge Capital at time t

	public static int[] theta = new int[50];				// int array of Compulsary level of education values chosen by government over generation 
	public static int[] theta_t = new int[1000];

	public static long[] LF_t = new long[1000];

// STRUCTURING THE " TECHNOLOGICAL CHANGE "

	public static float[] u_t = new float[1000];				// available technology state variable of the economy
	public static float[] udot_t = new float[1000];			// radical change in technology
	public static float[] u_upper_t = new float[1000];		// upgraded/modern technology state variable of the economy
	public static float[] u_lower_t = new float[1000];		// traditional technology state variable of the economy
	// u > 0; u_upper > u_lower; u_lower != 0;
	public static float[] c_ut_t = new float[1000];			// per capita set up cost of technology, function of both u(t) and t
	public static float[] c_t = new float[1000];				// per capita set up cost purely as a function of time

// THE ROLE OF THE GOVERNMENT

	public static float[] i_t = new float[1000];				// intermediate good
	public static float[] rd_t = new float[1000];			// resources devoted to RnD sector

	public static float[] h_t = new float[1000];				// quality adjusted human capital of an individual (human capital of a representative worker)
	public static float[] a_t = new float[1000];				// absorption capacity of the cohort
	public static float[] H_t = new float[1000];				// cohort specific human capital

	public static float[] G_1_t = new float[1000];			// total human capital of the first generation
	public static float[] G_2_t = new float[1000];			// total human capital of the second generation
	public static float[] G_3_t = new float[1000];			// total human capital of the third generation
	public static float[] v_t = new float[1000];			//present discounted value if all future incomes

	public static void main(String[] args){

	// FOUNDATIONAL THEORY MODEL 

		float add;
		K_t[0] = K_0;
		for (int i = 0; i < 4*T; i++){
			//generating KN_t values for -80 to 0 yrs to have a database to generate KN_t upon for future
			add = (float)Math.random()*(sum_var2)/10;
			sum_var2 = sum_var2 + add;
			KN_t[i] = sum_var2;

			if (i==0){
				Y_t[i] =Y_0;
			}
			else{
				//  Cobb Douglas type Production Function
				Y_t[i] = (float)Math.pow( K_t[i-1] , alpha ) * (float)Math.pow( KN_t[i-1], beta);		//  (1)
			}			

			// equation (3) discrete integration 
			K_t[i+1] = K_t[i] + s*Y_t[i] - (n+delta)*K_t[i];
			
		}		

	// EVOLUTIONARY VARIABLES
		B_t[0] = B_0;					// starting cohort ( also starting population )
		L_t[0] = B_t[0];				// starting population
		for (int i = 0; i < 4*T; i++){
			B_t[i+1] = (long)(B_t[i]*Math.exp(n*1));
			L_t[i+1] = L_t[i] + B_t[i+1];
		}

	//  Per Capita Cobb Douglas Production Function
		for (int i = 0; i < 4*T; i ++ ){
			k_t[i] = K_t[i]/L_t[i];
			kn_t[i] = KN_t[i]/L_t[i];
			y_t[i] = (float)Math.pow( k_t[i] , alpha ) *(float)Math.pow( kn_t[i], beta);			//  (2)	
		}


		int edu =1;
		for ( int i = 0; i < 4*T; i++){
			if (i%T == 0){
				// education value is set once every 20 years(generation lenght) and not every year
				edu = (int)(Math.random()*6)  + 6;
				theta[i/T] = edu;			
			}
			theta_t[i] = edu;							
		}
		

		long decay;
		int i=0;
		while (i <= theta[0]){
				LF_t[i] = LF_0;				//assumed Labour Force at start of base data
				i = i+1;
		}

		int k;
		for (i = 0; i < 4*T; i++){
			if (i < 3*T){
				decay = 0;
			}
			else{	
				decay = B_t[i%(3*T)];
			}

			if (LF_t[i+theta_t[i]-1] == 0){
				for (k=2;LF_t[i+theta_t[i]-k] == 0; k++){
				}
				long temp1 = LF_t[i+theta_t[i]-k];
				while(k>0){
					LF_t[i+theta_t[i]-k+1] = temp1;
					k--;
				}
				LF_t[i+theta_t[i]] = LF_t[i+theta_t[i]-1] + B_t[i] - decay;
			}
			else if (LF_t[i+theta_t[i]] != 0) {
				LF_t[i+theta_t[i]] = LF_t[i+theta_t[i]]  +B_t[i];
			}
			else{
				LF_t[i+theta_t[i]] = LF_t[i+theta_t[i]-1] + B_t[i] - decay;
			}	
		}

// STRUCTURING THE " TECHNOLOGICAL CHANGE "

		udot_t[0] = udot_0;						// size of the radical change in technology
		u_upper_t[0] = u_upper_0;
		u_lower_t[0] = u_lower_0;
		u_t[0] = u_0;		

		for (i = 1; i < 4*T; i++){
	// *** MAIN VARIABLE OF MODEL ***
			u_upper_t[i] = (float)Math.random()*3.0f + 2.0f + u_upper_t[i-1];
			u_lower_t[i] = (float)Math.random()*2.0f + u_lower_t[i-1];
			u_t[i] = u_upper_t[i] + u_lower_t[i];
			udot_t[i] = u_t[i] - u_t[i-1];
			c_t[i] = c_t[i-1] + ((float)Math.random()*0.1f)*c_t[i-1];
			// currently I am using random function to generate values but 
			// this will be needed to be taken care of 
			// Actually optimisation function for this is needed to be solved !!
		}
		
		for (i = 0; i < 4*T; i++){
			if ( udot_t[i] > 0){
				c_ut_t[i] = c_t[i] * (u_upper_t[i] - u_lower_t[i]);		//function of both u(t) and t
			}
			else{
				c_ut_t[i] = 0.0f;			// set up cost only in case of radical change
			}
		}

	// THE ROLE OF THE GOVERNMENT
		
		for (i = 0; i < 4*T; i++){
			rd_t[i] = eta * y_t[i];
			i_t[i] = u_t[i]*rd_t[i];									//  (7)

			//i_t = eta * u_t * y_t;									//  (8)
			//the cohort specific human capital then becomes
			a_t[i] = phi * theta_t[i];
			h_t[i] = a_t[i] * i_t[i];									//   (9)
			H_t[i] = h_t[i] * B_t[i];
			
		}


		int t = 4*T;				// primary time variable
		int j; 

		float sum1, sum2, sum3, sum4;
		float psi = psi_const;							//human capital combination exponent for first generation
		float gamma = gamma_const;						//human capital combination exponent for second generation

		for (j=0;j<20*T;j++){
			B_t[t+j+1] = (long)(B_t[t+j] * Math.exp(n*1));
			L_t[t+j+1] = L_t[t+j] + B_t[t+j] - B_t[j];
		}

		for(j=-T;j<=0;j++){
			sum1=0.0f;
			sum2 = 0.0f; 
			sum3=0.0f;
			for (i=0; i<T;i++){
				if (i >= theta_t[t+j]){
					sum1 = sum1 + H_t[t+j+theta_t[t+j]-i] * B_t[t+j + theta_t[t+j]-i];
				}
				sum2 = sum2 + H_t[t+j+theta_t[t+j]-i-T] * B_t[t+j + theta_t[t+j] -i-T];
				sum3 = sum3 + H_t[t+j+theta_t[t+j]-i-2*T] * B_t[t+j + theta_t[t+j] -i-2*T];
			} 
			G_1_t[t+j+ theta_t[t+j]] = sum1;
			G_2_t[t+j+ theta_t[t+j]] = sum2;
			G_3_t[t+j+ theta_t[t+j]] = sum3;
			// the firms combine the human capital of all three generations of indivaduals in the labour force
			// to produce an economy wide knowledge capital which is input for aggregate production function
			KN_t[t+j+ theta_t[t+j]] = (float)Math.pow(G_1_t[t+j+ theta_t[t+j]], psi) 
								*(float) Math.pow(G_2_t[t+j+ theta_t[t+j]], gamma) 
											* (float)Math.pow( G_3_t[t+j+ theta_t[t+j]], 1 - psi - gamma);


			Y_t[t+j+theta_t[t+j]] = (float)Math.pow( K_t[t+j+theta_t[t+j]-1] , alpha ) * (float)Math.pow( KN_t[t+j+theta_t[t+j]-1], beta);		//  (1)
			K_t[t+j+theta_t[t+j]+1] = K_t[t+j+theta_t[t+j]] + s*Y_t[t+j+theta_t[t+j]] - (n+delta)*K_t[t+j+theta_t[t+j]];

			k_t[t+j+theta_t[t+j]+1] = K_t[t+j+theta_t[t+j]+1]/L_t[t+j+theta_t[t+j]+1];
			kn_t[t+j+theta_t[t+j]] = KN_t[t+j+theta_t[t+j]]/L_t[t+j+theta_t[t+j]];
			y_t[t+j+theta_t[t+j]] = (float)Math.pow( k_t[t+j+theta_t[t+j]-1] , alpha ) *(float)Math.pow( kn_t[t+j+theta_t[t+j]-1], beta);

		}

// ALL Evolutionary BASE DATA STORED ARE APPROPRIATELY GENERATED BY NOW FROM T = -80 TO T = -1 YRS
// NOW T = 0 YRS ONWARDS WE"LL GENERATE (here T = 0 means array index 80 and onwards)


		for(j=0;j<20*T;j++){

			if ((t+j)%T == 0){
				// education value is set once every 20 years(generation lenght) and not every year
				edu = (int)(Math.random()*6)  + 6;
	// *** MAIN VARIABLE OF MODEL ***
				theta[(t+j)/T] = edu;			
			}
			theta_t[t+j] = edu;

			if (j<350){
				sum4 = 0.0f;
				for (k=0;k<theta_t[t+j];k++){
					sum4 = sum4 + (float)(y_t[t+j+k]*Math.exp(-rho*k));
				}

				for (k=theta_t[t+j];k<100;k++){
					sum4 = sum4+(float)(y_t[t+j+theta_t[t+j]-1]*Math.exp(-rho*k));
				}
				//System.out.println(y_t[t+j+theta_t[t+j]-1]);
				//System.out.println(sum4 + "!!");
				v_t[t+j] = sum4;
			}


			if (LF_t[t+j+theta_t[t+j]-1] == 0){
				// skipped year i.e the expected labour force needs to do more years because
				// of increased cumpulsary level of education
				for (i=2;LF_t[t+j+theta_t[t+j]-i] == 0; i++){
				}
				long temp1 = LF_t[t+j+theta_t[t+j]-i];
				float temp2 = KN_t[t+j+ theta_t[t+j]-i];
				int temp3;
				while(i>0){
					LF_t[t+j+theta_t[t+j]-i+1] = temp1;
					// for the interim years the labour force remains constant
					// and so does the knowledge capital
					KN_t[t+j+ theta_t[t+j]-i+1] = temp2;
					Y_t[t+j+theta_t[t+j]-i+1] = (float)Math.pow( K_t[t+j+theta_t[t+j]-i] , alpha ) * (float)Math.pow( KN_t[t+j+theta_t[t+j]-i], beta);		//  (1)
					K_t[t+j+theta_t[t+j]-i+1] = K_t[t+j+theta_t[t+j]-i] + s*Y_t[t+j+theta_t[t+j]-i] - (n+delta)*K_t[t+j+theta_t[t+j]-i];
					k_t[t+j+theta_t[t+j]-i+1] = K_t[t+j+theta_t[t+j]-i+1]/L_t[t+j+theta_t[t+j]-i+1];
					kn_t[t+j+theta_t[t+j]-i+1] = KN_t[t+j+theta_t[t+j]-i+1]/L_t[t+j+theta_t[t+j]-i];
					y_t[t+j+theta_t[t+j]-i] = (float)Math.pow( k_t[t+j+theta_t[t+j]-i-1] , alpha ) *(float)Math.pow( kn_t[t+j+theta_t[t+j]-i-1], beta);

					i--;
				}
				LF_t[t+j+theta_t[t+j]] = LF_t[t+j+theta_t[t+j]-1] + B_t[t+j] - B_t[(t+j)-(3*T)];

				sum1=0.0f;
				sum2 = 0.0f; 
				sum3=0.0f;
				for (i=0; i<T;i++){
					if (i >= theta_t[t+j]){
						sum1 = sum1 + H_t[t+j+theta_t[t+j]-i] * B_t[t+j + theta_t[t+j]-i];
					}
					sum2 = sum2 + H_t[t+j+theta_t[t+j]-i-T] * B_t[t+j + theta_t[t+j] -i-T];
					sum3 = sum3 + H_t[t+j+theta_t[t+j]-i-2*T] * B_t[t+j + theta_t[t+j] -i-2*T];
				} 
				G_1_t[t+j+ theta_t[t+j]] = sum1;
				G_2_t[t+j+ theta_t[t+j]] = sum2;
				G_3_t[t+j+ theta_t[t+j]] = sum3;
				// System.out.println(G_1_t[index] + ",G1(t): " + index);
				// the firms combine the human capital of all three generations of indivaduals in the labour force
				// to produce an economy wide knowledge capital which is input for aggregate production function
				KN_t[t+j+ theta_t[t+j]] = (float)(Math.pow(G_1_t[t+j+ theta_t[t+j]], psi) 
									* Math.pow(G_2_t[t+j+ theta_t[t+j]], gamma) 
												* Math.pow( G_3_t[t+j+ theta_t[t+j]], 1 - psi - gamma));

			}
			else if (LF_t[t+j+theta_t[t+j]] != 0) {
				//crowded year with more than 1 cohort joining the labour force due to decreased level of 
				//cumpulsary education for future cohorts
				LF_t[t+j+theta_t[t+j]] = LF_t[t+j+theta_t[t+j]] +B_t[t+j];
			}
			else{
				LF_t[t+j+theta_t[t+j]] = LF_t[t+j+theta_t[t+j]-1] + B_t[t+j] - B_t[(t+j)-(3*T)];

				sum1=0.0f;
				sum2 = 0.0f; 
				sum3=0.0f;
				for (i=0; i<T;i++){
					if (i >= theta_t[t+j]){
						sum1 = sum1 + H_t[t+j+theta_t[t+j]-i] * B_t[t+j + theta_t[t+j]-i];
					}
					sum2 = sum2 + H_t[t+j+theta_t[t+j]-i-T] * B_t[t+j + theta_t[t+j] -i-T];
					sum3 = sum3 + H_t[t+j+theta_t[t+j]-i-2*T] * B_t[t+j + theta_t[t+j] -i-2*T];
				} 
				G_1_t[t+j+ theta_t[t+j]] = sum1;
				G_2_t[t+j+ theta_t[t+j]] = sum2;
				G_3_t[t+j+ theta_t[t+j]] = sum3;
				// System.out.println(G_1_t[index] + ",G1(t): " + index);
				// the firms combine the human capital of all three generations of indivaduals in the labour force
				// to produce an economy wide knowledge capital which is input for aggregate production function
				KN_t[t+j+ theta_t[t+j]] = (float)(Math.pow(G_1_t[t+j+ theta_t[t+j]], psi) 
									* Math.pow(G_2_t[t+j+ theta_t[t+j]], gamma) 
												* Math.pow( G_3_t[t+j+ theta_t[t+j]], 1 - psi - gamma));

				Y_t[t+j+theta_t[t+j]] = (float)Math.pow( K_t[t+j+theta_t[t+j]-1] , alpha ) * (float)Math.pow( KN_t[t+j+theta_t[t+j]-1], beta);		//  (1)
				K_t[t+j+theta_t[t+j]+1] = K_t[t+j+theta_t[t+j]] + s*Y_t[t+j+theta_t[t+j]] - (n+delta)*K_t[t+j+theta_t[t+j]];

				// System.out.println(K_t[index] + ",K_t(t+j): " + index);
				// System.out.println(KN_t[index] + ",KN_t(t+j): " + index);
				// System.out.println(LF_t[index] + ",LF(t+j): " + index);
				

				k_t[t+j+theta_t[t+j]+1] = K_t[t+j+theta_t[t+j]+1]/L_t[t+j+theta_t[t+j+1]];
				kn_t[t+j+theta_t[t+j]] = KN_t[t+j+theta_t[t+j]]/L_t[t+j+theta_t[t+j]];
				y_t[t+j+theta_t[t+j]] = (float)Math.pow( k_t[t+j+theta_t[t+j]-1] , alpha ) *(float)Math.pow( kn_t[t+j+theta_t[t+j]-1], beta);

			}


			
	// *** MAIN VARIABLE OF MODEL ***
			u_upper_t[t+j] = (float)Math.random()*25 + 15.0f + u_upper_t[t+j-1];
			u_lower_t[t+j] = (float)Math.random()*15 + u_lower_t[t+j-1];
			u_t[t+j] = u_upper_t[t+j] + u_lower_t[t+j];
			udot_t[t+j] = u_t[t+j] - u_t[t+j-1];
			c_t[t+j] = (float)Math.random()*10 + 10.0f;
			// currently I am using random function to generate values but 
			// this will be needed to be taken care of 
			// Actually optimisation function for this is needed to be solved !!
		
			if ( udot_t[t+j] > 0){
				c_ut_t[t+j] = c_t[t+j] * (u_upper_t[t+j] - u_lower_t[t+j]);		//function of both u(t) and t
			}
			else{
				c_ut_t[t+j] = 0.0f;			// set up cost only in case of radical change
			}

			rd_t[t+j] = eta * y_t[t+j];
			i_t[t+j] = u_t[t+j]*rd_t[t+j];									//  (7)

			//i_t = eta * u_t * y_t;									//  (8)
			//the cohort specific human capital then becomes
			a_t[t+j] = phi * theta_t[t+j];
			h_t[t+j] = a_t[t+j] * i_t[t+j];									//   (9)
			H_t[t+j] = h_t[t+j] * B_t[t+j];
			
		}	


// ==========|||||||||||||||||====================||||||||||||===================|||||||||||||||============================|||||||||||

		// //Above here we generate a random timeseries data array for physical annual capital
		// System.out.println("1) K(t): " + Arrays.toString(K_t)); System.out.println("\n\n ");
		// //Above here we generate a random timeseries data array for knowledge annual capital
		// System.out.println("2) KN(t): " + Arrays.toString(KN_t)); System.out.println("\n\n ");
		// //Above here we calculate timeseries data array for GDP based on generated values of K(t) and KN(t)
		// System.out.println("3) Y(t): " + Arrays.toString(Y_t)); System.out.println("\n\n ");

		// //Above here we calculate timeseries data array for cohort sizes based on first cohort size
		// System.out.println("4) B(t): " + Arrays.toString(B_t)); System.out.println("\n\n ");
		// //Above here we generate timeseries data array for population for first 80 yrs to act as our base population data
		// System.out.println("5) L(t): " + Arrays.toString(L_t)); System.out.println("\n\n ");
		// //Above here we compute timeseries data array for Labour Force
		// System.out.println("6) LF(t): " + Arrays.toString(LF_t)); System.out.println("\n\n ");

		// //Above here we generate per capita timeseries data array for k(t), kn(t), y(t)
		// System.out.println("7) k(t): " + Arrays.toString(k_t)); System.out.println("\n\n ");
		// System.out.println("8) kn(t): " + Arrays.toString(kn_t)); System.out.println("\n\n ");
		// System.out.println("9) y(t): " + Arrays.toString(y_t)); System.out.println("\n\n ");

		// //Above here we generate timeseries data array for compulsary education level
		// System.out.println("10)a) theta: " + Arrays.toString(theta)); System.out.println("\n\n ");
		// System.out.println("10)b) theta(t): " + Arrays.toString(theta_t)); System.out.println("\n\n ");
		
		// //Above here we start computing timeseries data array for technology variables
		// System.out.println("11) u_upper(t): " + Arrays.toString(u_upper_t)); System.out.println("\n\n ");
		// System.out.println("12) u_lower(t): " + Arrays.toString(u_lower_t)); System.out.println("\n\n ");
		// System.out.println("13) u(t): " + Arrays.toString(u_t)); System.out.println("\n\n ");
		// System.out.println("14) udot_t(t): "  + Arrays.toString(udot_t)); System.out.println("\n\n ");
		// System.out.println("15) c(t): "  + Arrays.toString(c_t)); System.out.println("\n\n ");

		// //Above here we start computing timeseries data array for set up cost
		// System.out.println("16) c_u(t): " + Arrays.toString(c_ut_t)); System.out.println("\n\n ");

		// //Above here we start computing timeseries data array for income and RnD investment
		// System.out.println("17) rd(t): " + Arrays.toString(rd_t)); System.out.println("\n\n ");
		// System.out.println("18) i(t): " + Arrays.toString(i_t)); System.out.println("\n\n ");

		// //Above here we start computing timeseries data array for knowledge variables
		// System.out.println("19) a(t): " + Arrays.toString(a_t)); System.out.println("\n\n ");
		// System.out.println("20) h(t): " + Arrays.toString(h_t)); System.out.println("\n\n ");
		// System.out.println("21) H(t): " + Arrays.toString(H_t)); System.out.println("\n\n ");
		
		// //Above here we start computing timeseries data array for Cohort Sizes
		// System.out.println("22) B(t): " + Arrays.toString(B_t)); System.out.println("\n\n ");
		// //Above here we start computing timeseries data array for Population
		// System.out.println("23) L(t): " + Arrays.toString(L_t)); System.out.println("\n\n ");
		// //Above here we start computing timeseries data array for Population
		// System.out.println("24) LF(t): " + Arrays.toString(LF_t)); System.out.println("\n\n ");		

		// System.out.println("25) G_1(t): " + Arrays.toString(G_1_t)); System.out.println("\n\n ");
		// System.out.println("26) G_2(t): " + Arrays.toString(G_2_t)); System.out.println("\n\n ");
		// System.out.println("27) G_3(t): " + Arrays.toString(G_3_t)); System.out.println("\n\n ");
		
		// System.out.println("28) KN_t(t): " + Arrays.toString(KN_t)); System.out.println("\n\n ");
		

// ==========|||||||||||||||||====================||||||||||||===================|||||||||||||||============================|||||||||||


		XYLineChart_AWT_float chart1 = new XYLineChart_AWT_float(" Figure", "Y(t) = K(t)^alpha * KN(t)^beta", "Y(t)", Y_t, 20);
		chart1.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart1 );          
		chart1.setVisible( true );

		XYLineChart_AWT_float chart2 = new XYLineChart_AWT_float(" Figure", "K(t)", "K(t)", K_t, 20);
		chart2.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart2 );          
		chart2.setVisible( true );

		XYLineChart_AWT_float chart3 = new XYLineChart_AWT_float(" Figure", "KN(t)", "KN(t)", KN_t, 20);
		chart3.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart3 );          
		chart3.setVisible( true );

		XYLineChart_AWT_long chart4 = new XYLineChart_AWT_long(" Figure", "B(t)", "B(t)", B_t, 20);
		chart4.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart4 );          
		chart4.setVisible( true );

		XYLineChart_AWT_long chart5 = new XYLineChart_AWT_long(" Figure", "L(t)", "L(t)", L_t, 20);
		chart5.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart5 );          
		chart5.setVisible( true );

		XYLineChart_AWT_long chart6 = new XYLineChart_AWT_long(" Figure", "LF(t)", "LF(t)", LF_t, 20);
		chart6.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart6 );          
		chart6.setVisible( true );

		XYLineChart_AWT_float chart7 = new XYLineChart_AWT_float(" Figure", "k(t)", "k(t)", k_t, 20);
		chart7.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart7 );          
		chart7.setVisible( true );

		XYLineChart_AWT_float chart8 = new XYLineChart_AWT_float(" Figure", "y(t)", "y(t)", y_t, 20);
		chart8.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart8 );          
		chart8.setVisible( true );

		XYLineChart_AWT_float chart9 = new XYLineChart_AWT_float(" Figure", "kn(t)", "kn(t)", kn_t, 20);
		chart9.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart9 );          
		chart9.setVisible( true );

		XYLineChart_AWT_int chart10 = new XYLineChart_AWT_int(" Figure", "theta_t(t)", "theta_t(t)", theta_t, 20);
		chart10.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart10 );          
		chart10.setVisible( true );

		XYLineChart_AWT_float chart11 = new XYLineChart_AWT_float(" Figure", "u_upper(t)", "u_upper(t)", u_upper_t, 20);
		chart11.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart11 );          
		chart11.setVisible( true );

		XYLineChart_AWT_float chart12 = new XYLineChart_AWT_float(" Figure", "u_lower(t)", "u_lower(t)", u_lower_t, 20);
		chart12.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart12 );          
		chart12.setVisible( true );

		XYLineChart_AWT_float chart13 = new XYLineChart_AWT_float(" Figure", "u(t)", "u(t)", u_t, 20);
		chart13.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart13 );          
		chart13.setVisible( true );

		XYLineChart_AWT_float chart14 = new XYLineChart_AWT_float(" Figure", "c(t)", "c(t)", c_t, 20);
		chart14.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart14 );          
		chart14.setVisible( true );

		XYLineChart_AWT_float chart15 = new XYLineChart_AWT_float(" Figure", "c_ut(t)", "c_ut(t)", c_ut_t, 20);
		chart15.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart15 );          
		chart15.setVisible( true );

		XYLineChart_AWT_float chart16 = new XYLineChart_AWT_float(" Figure", "i(t)", "i(t)", i_t, 20);
		chart16.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart16 );          
		chart16.setVisible( true );

		XYLineChart_AWT_float chart17 = new XYLineChart_AWT_float(" Figure", "a(t)", "a(t)", a_t, 20);
		chart17.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart17 );          
		chart17.setVisible( true );

		XYLineChart_AWT_float chart18 = new XYLineChart_AWT_float(" Figure", "h(t)", "h(t)", h_t, 20);
		chart18.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart18 );          
		chart18.setVisible( true );

		XYLineChart_AWT_float chart19 = new XYLineChart_AWT_float(" Figure", "rd(t)", "rd(t)", rd_t, 20);
		chart19.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart19 );          
		chart19.setVisible( true );

		XYLineChart_AWT_float chart20 = new XYLineChart_AWT_float(" Figure", "H(t)", "H(t)", H_t, 20);
		chart20.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart20 );          
		chart20.setVisible( true );

		XYLineChart_AWT_float chart21 = new XYLineChart_AWT_float(" Figure", "v(t)", "v(t)", v_t, 20);
		chart21.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart21 );          
		chart21.setVisible( true );
	}
	
}