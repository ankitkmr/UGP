import java.awt.Color; 
import java.awt.BasicStroke; 
import org.jfree.chart.ChartPanel; 
import org.jfree.chart.JFreeChart; 
import org.jfree.data.xy.XYDataset; 
import org.jfree.data.xy.XYSeries; 
import org.jfree.ui.ApplicationFrame; 
import org.jfree.ui.RefineryUtilities; 
import org.jfree.chart.plot.XYPlot; 
import org.jfree.chart.ChartFactory; 
import org.jfree.chart.plot.PlotOrientation; 
import org.jfree.data.xy.XYSeriesCollection; 
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;

import java.util.Arrays;


public class TheoryModel{

		////  ** Arbitrary Value Taken Here ** ////

		// FOUNDATIONAL THEORY MODEL 
		static int T = 20 ;			// number of years in a generation

		static float K_0 = 100000.0f;				// starting physical capital
		static float sum_var2 = 30000.0f ;			// starting knowledge capital

		static float alpha = 0.6f;
		static float beta = 0.4f;		//  0 < alpha, beta < 1

		static float s = 0.3f;					// constant fraction of output that is saved and invested
		static float delta = 0.1f;				// exogenously fixed capital depreciation rate
		static float n = 0.18f;				// exogenous growth rate of cohort
	
		// THE ROLE OF THE GOVERNMENT
		static float eta = 0.4f;					// constant fraction of income y\(t\) that is invested in the RnD
		static float phi = 0.35f;					// phi > 0, conversion rate of education to cohort specific absorpition capacity 

		static float udot_0 = 0.3f;						// size of the radical change in technology
		static float u_upper_0 = 60.0f;
		static float u_lower_0 = 40.0f;
		static float u_0 = 100.0f;	
		////  ** Arbitrary Value Taken Here ** ////

	// FOUNDATIONAL THEORY MODEL 

		//  Cobb Douglas type Production Function variables
		public static float[] Y_t = new float[300];		// GDP of Economy at time t
		public static float[] K_t = new float[300];		// Physical Capital at time t
		public static float[] KN_t = new float[300];		// Knowledge Capital at time t

	// EVOLUTIONARY VARIABLES

		public static long B_0 = 1000;
		public static long[] B_t = new long[300];		// new Cohort born at each time t ( number of people born at time t )
		public static long[] L_t = new long[300];		// size of the population at time 	

		public static float[] y_t = new float[300];		// Per Capita GDP of Economy at time t
		public static float[] k_t = new float[300];		// Physical Capital per worker at time t
		public static float[] kn_t = new float[300];		// Per Capita Knowledge Capital at time t
		public static float[] k_dot_t = new float[300];	// evolution rate of physical capital overtime

		public static int[] theta = new int[15];			// int array of Compulsary level of education values chosen by government over generation 
		public static int[] theta_t = new int[300];

		public static long[] LF_t = new long[300];

	// STRUCTURING THE " TECHNOLOGICAL CHANGE "

		public static float[] u_t = new float[300];			// available technology state variable of the economy
		public static float[] udot_t = new float[300];		// radical change in technology
		public static float[] u_upper_t = new float[300];		// upgraded/modern technology state variable of the economy
		public static float[] u_lower_t = new float[300];		// traditional technology state variable of the economy
		// u > 0; u_upper > u_lower; u_lower != 0;
		public static float[] c_ut_t = new float[300];			// per capita set up cost of technology, function of both u(t) and t
		public static float[] c_t = new float[300];				// per capita set up cost purely as a function of time

		static int tau = 80;					// secondary time variable

	// THE ROLE OF THE GOVERNMENT

		public static float[] i_t = new float[300];		// intermediate good
		public static float[] rd_t = new float[300];		// resources devoted to RnD sector

		public static float[] h_t = new float[300];		// quality adjusted human capital of an individual (human capital of a representative worker)
		public static float[] a_t = new float[300];		// absorption capacity of the cohort
		public static float[] H_t = new float[300];		// cohort specific human capital



	public static void main(String[] args){

// FOUNDATIONAL THEORY MODEL 

		float add;
		for (int i = 0; i < 4*T; i++){
			//generating KN_t values for -80 to 0 yrs to have a database to generate KN_t upon for future
			add = (float)Math.random()*100.0f + 1;
			sum_var2 = sum_var2 + add;
			KN_t[i] = sum_var2;

			if (i == 0){
				K_t[i] = K_0;
			}
			else{
				// equation (3) discrete integration 
				K_t[i] = K_t[i-1] + s*Y[i-1] - (n+delta)*K_t[i-1];
			}

			//  Cobb Douglas type Production Function
			Y_t[i] = (float)Math.pow( K_t[i] , alpha ) * (float)Math.pow( KN_t[i], beta);		//  (1)
		}		


// EVOLUTIONARY VARIABLES
		B_t[0] = B_0;					// starting cohort ( also starting population )
		L_t[0] = B_t[0];				// starting population
		for (int i = 1; i < 4*T; i++){
			B_t[i] = (long)(B_t[i-1]*Math.exp(n*1));
			L_t[i] = L_t[i-1] + B_t[i];
		}


	//  Per Capita Cobb Douglas Production Function
		for (int i = 0; i < 4*T; i ++ ){
			k_t[i] = K_t[i]/L_t[i];
			kn_t[i] = KN_t[i]/L_t[i];
			y_t[i] = (float)Math.pow( k_t[i] , alpha ) *(float)Math.pow( kn_t[i], beta);			//  (2)	
		}

		//Above here we generate a random timeseries data array for physical annual capital
		System.out.println("1) K(t): " + Arrays.toString(K_t)); System.out.println("\n\n ");
		//Above here we generate a random timeseries data array for knowledge annual capital
		System.out.println("2) KN(t): " + Arrays.toString(KN_t)); System.out.println("\n\n ");
		//Above here we calculate timeseries data array for GDP based on generated values of K(t) and KN(t)
		System.out.println("3) Y(t): " + Arrays.toString(Y_t)); System.out.println("\n\n ");

		//Above here we calculate timeseries data array for cohort sizes based on first cohort size
		System.out.println("4) B(t): " + Arrays.toString(B_t)); System.out.println("\n\n ");
		//Above here we generate timeseries data array for population for first 80 yrs to act as our base population data
		System.out.println("5) L(t): " + Arrays.toString(L_t)); System.out.println("\n\n ");

		//Above here we generate per capita timeseries data array for k(t), kn(t), y(t)
		System.out.println("6) k(t): " + Arrays.toString(k_t)); System.out.println("\n\n ");
		System.out.println("7) kn(t): " + Arrays.toString(kn_t)); System.out.println("\n\n ");
		System.out.println("8) y(t): " + Arrays.toString(y_t)); System.out.println("\n\n ");


		// //  evolution of physical capital 
		// for (int i = 0; i <4*T; i++){
		// 	k_dot_t[i] = s * y_t[i] - (n + delta) * k_t[i];							//  (3)
		// }
		// //Above here we compute timeseries data array for evolution physical capital
		// System.out.println("9) k_dot(t): " + Arrays.toString(k_dot_t)); System.out.println("\n\n ");

		int edu = 1;
		for ( int i = 0; i < 4*T; i++){
			if (i == 0 || i == T || i == 2*T || i == 3*T){
				// education value is set once every 20 years(generation lenght) and not every year
				edu = (int)(Math.random()*6)  + 6;
	// *** MAIN VARIABLE OF MODEL ***
				theta[(i+1)/T] = edu;			
			}

			theta_t[i] = edu;							
		}
		//Above here we generate timeseries data array for compulsary education level
		System.out.println("10)a) theta: " + Arrays.toString(theta)); System.out.println("\n\n ");
		//Above here we generate timeseries data array for compulsary education level
		System.out.println("10)b) theta(t): " + Arrays.toString(theta_t)); System.out.println("\n\n ");

		long decay;
		int i=0;
		while (i <= theta[0]){
				LF_t[i] = 10;				//10 assumed Labour Force at start of base data
				i = i+1;
		}

		for (i = 0; i < 4*T; i++){
			if (i < 3*T){
				decay = 0;
			}
			else{	
				decay = B_t[(i+1)%(3*T)];
			}
			LF_t[i+theta[(i+1)/T]+1] = LF_t[i+theta[(i+1)/T]] + B_t[i] - decay;	
		}

		//Above here we compute timeseries data array for Labour Force
		System.out.println("11) LF(t): " + Arrays.toString(LF_t)); System.out.println("\n\n ");


// ALL Evolutionary BASE DATA STORED AND APPROPRIATELY GENERATED BY NOW 
// FROM T = -80 TO T = -1 YRS

// NOW T = 0 YRS ONWARDS WE"LL GENERATE (here T = 0 means array index 80 and onwards)

		int t = 4*T - 1;				// primary time variable

		// Cohort size function
		// Generating cohort sizes from T = t + 1 to T = tau 
		for (i = 0; i < tau; i++){
			B_t[t+i] = (long)(B_t[t+i-1] * Math.exp(n));								//  (4)
		}
		//Above here we start computing timeseries data array for Cohort Sizes
		System.out.println("12) B(t): " + Arrays.toString(B_t)); System.out.println("\n\n ");

		//  Total Population at any time
		//  Generating Population sizes from T = t + 1 to T = tau 
		for (i = 0; i < tau; i++){
			L_t[t+i] = (long)(L_t[t+i-1] + B_t[t+i] - B_t[i]);
		}
		//Above here we start computing timeseries data array for Population
		System.out.println("13) L(t): " + Arrays.toString(L_t)); System.out.println("\n\n ");

		//Total Labour Force Participation at any time t
		for (i = 0; i < 4 * T; i++){
			LF_t[t+i+theta_t[(i)/T]] = LF_t[t+i+theta_t[(i)/T]-1] + B_t[t+i] - B_t[(i)%(3*T)];
		}
		//Above here we start computing timeseries data array for Population
		System.out.println("14) LF(t): " + Arrays.toString(LF_t)); System.out.println("\n\n ");


	// STRUCTURING THE " TECHNOLOGICAL CHANGE "

		udot_t[0] = udot_0;						// size of the radical change in technology
		u_upper_t[0] = u_upper_0;
		u_lower_t[0] = u_lower_0;
		u_t[0] = u_0;		

		for (i = 1; i < 4*T; i++){
// *** MAIN VARIABLE OF MODEL ***
			u_upper_t[i] = (float)Math.random()*25 + 15.0f + u_upper_t[i-1];
			u_lower_t[i] = (float)Math.random()*15 + u_lower_t[i-1];
			u_t[i] = u_upper_t[i] + u_lower_t[i];
			udot_t[i] = u_t[i] - u_t[i-1];
			c_t[i] = (float)Math.random()*10 + 10.0f;
			// currently I am using random function to generate values but 
			// this will be needed to be taken care of 
			// Actually optimisation function for this is needed to be solved !!
		
		}
		//Above here we start computing timeseries data array for technology variables
		System.out.println("15) u_upper(t): " + Arrays.toString(u_upper_t)); System.out.println("\n\n ");
		System.out.println("16) u_lower(t): " + Arrays.toString(u_lower_t)); System.out.println("\n\n ");
		System.out.println("17) u(t): " + Arrays.toString(u_t)); System.out.println("\n\n ");
		System.out.println("18) udot_t(t): "  + Arrays.toString(udot_t)); System.out.println("\n\n ");
		System.out.println("19) c(t): "  + Arrays.toString(c_t)); System.out.println("\n\n ");

		for (i = 0; i < 4*T; i++){
			if ( udot_t[i] > 0){
				c_ut_t[i] = c_t[i] * (u_upper_t[i] - u_lower_t[i]);		//function of both u(t) and t
			}
			else{
				c_ut_t[i] = 0.0f;			// set up cost only in case of radical change
			}
		}
		//Above here we start computing timeseries data array for set up cost
		System.out.println("20) c_u(t): " + Arrays.toString(c_ut_t)); System.out.println("\n\n ");

	// THE ROLE OF THE GOVERNMENT
		
		for (i = 0; i < 4*T; i++){
			rd_t[i] = eta * y_t[i];
			i_t[i] = u_t[i]*rd_t[i];									//  (7)
			
		}
		//Above here we start computing timeseries data array for income and RnD investment
		System.out.println("21) rd(t): " + Arrays.toString(rd_t)); System.out.println("\n\n ");
		System.out.println("22) i(t): " + Arrays.toString(i_t)); System.out.println("\n\n ");

		//i_t = eta * u_t * y_t;											//  (8)

		for (i = 0; i < 4*T; i++){
			//the cohort specific human capital then becomes
			a_t[i] = phi * theta_t[i];
			h_t[i] = a_t[i] * i_t[i];												//   (9)
			H_t[i] = h_t[i] * B_t[i];
		}
		//Above here we start computing timeseries data array for knowledge variables
		System.out.println("23) a(t): " + Arrays.toString(a_t)); System.out.println("\n\n ");
		System.out.println("24) h(t): " + Arrays.toString(h_t)); System.out.println("\n\n ");
		System.out.println("25) H(t): " + Arrays.toString(H_t)); System.out.println("\n\n ");
		
		// **** NOT Required **** //
		
		// float[] h_t_years;			// float array of h(t) values recorded overtime at the start of each year
		// float[] B_t_years;			// float array of B(t) values recorded overtime at the start of each year
		
		// **** NOT Required **** //

		float G_1_t_plus_theta = 10.0f;		// total human capital of the first generation

		// for ( float sum = 0, int i = 0; i <= T; i = i + 1){
		// 	sum = sum + h_t[t + theta[t] - i] * B_t[t + theta[t] - i] * 1;
		// }

		// G_1_t_plus_theta = sum

		float G_2_t_plus_theta = 10.0f;		// total human capital of the second generation

		// for ( sum = 0, i = T; i <= 2 * T; i = i + 1){
		// 	sum = sum + h_t[t + theta[t] - i] * B_t[t + theta[t] - i] * 1;
		// }

		// G_2_t_plus_theta = sum

		float G_3_t_plus_theta = 10.0f;		// total human capital of the third generation

		// for ( sum = 0, i = 2 * T; i <= 3 * T; i = i + 1){
		// 	sum = sum + h_t[t + theta[t] - i] * B_t[t + theta[t] - i] * 1;
		// }

		// G_3_t_plus_theta = sum

		float psi = 0.2f;			//human capital combination eiponent for first generation
		float gamma = 0.3f;		//human capital combination eiponent for second generation

		// the firms combine the human capital of all three generations of indivaduals in the labour force
		// to produce an economy wide knowledge capital which is input for aggregate production function
		float KN_t_plus_theta = (float)Math.pow(G_1_t_plus_theta, psi) 
										*(float) Math.pow(G_2_t_plus_theta, gamma) 
												* (float)Math.pow( G_3_t_plus_theta, 1 - psi - gamma); 

		XYLineChart_AWT_float chart1 = new XYLineChart_AWT_float(" Figure", "Y(t) = K(t)^alpha * KN(t)^beta", "Y(t)", Y_t, 4);
      	chart1.pack( );          
      	RefineryUtilities.centerFrameOnScreen( chart1 );          
      	chart1.setVisible( true );

      	XYLineChart_AWT_float chart2 = new XYLineChart_AWT_float(" Figure", "Y(t) = K(t)^alpha * KN(t)^beta", "K(t)", K_t, 4);
      	chart2.pack( );          
      	RefineryUtilities.centerFrameOnScreen( chart2 );          
      	chart2.setVisible( true );

      	XYLineChart_AWT_float chart3 = new XYLineChart_AWT_float(" Figure", "Y(t) = K(t)^alpha * KN(t)^beta", "KN(t)", KN_t, 4);
      	chart3.pack( );          
      	RefineryUtilities.centerFrameOnScreen( chart3 );          
      	chart3.setVisible( true );

		XYLineChart_AWT_long chart4 = new XYLineChart_AWT_long(" Figure", "B(t)", "B(t)", B_t, 8);
		chart4.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart4 );          
		chart4.setVisible( true );

		XYLineChart_AWT_long chart5 = new XYLineChart_AWT_long(" Figure", "L(t)", "L(t)", L_t, 8);
		chart5.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart5 );          
		chart5.setVisible( true );

		XYLineChart_AWT_long chart6 = new XYLineChart_AWT_long(" Figure", "LF(t)", "LF(t)", LF_t, 8);
		chart6.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart6 );          
		chart6.setVisible( true );

		XYLineChart_AWT_float chart7 = new XYLineChart_AWT_float(" Figure", "k(t)", "k(t)", k_t, 4);
		chart7.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart7 );          
		chart7.setVisible( true );

		XYLineChart_AWT_float chart8 = new XYLineChart_AWT_float(" Figure", "y(t)", "y(t)", y_t, 4);
		chart8.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart8 );          
		chart8.setVisible( true );

		XYLineChart_AWT_float chart9 = new XYLineChart_AWT_float(" Figure", "kn(t)", "kn(t)", kn_t, 4);
		chart9.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart9 );          
		chart9.setVisible( true );

		XYLineChart_AWT_int chart10 = new XYLineChart_AWT_int(" Figure", "theta_t(t)", "theta_t(t)", theta_t, 4);
		chart10.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart10 );          
		chart10.setVisible( true );

		XYLineChart_AWT_float chart11 = new XYLineChart_AWT_float(" Figure", "u_upper(t)", "u_upper(t)", u_upper_t, 4);
		chart11.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart11 );          
		chart11.setVisible( true );

		XYLineChart_AWT_float chart12 = new XYLineChart_AWT_float(" Figure", "u_lower(t)", "u_lower(t)", u_lower_t, 4);
		chart12.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart12 );          
		chart12.setVisible( true );

		XYLineChart_AWT_float chart13 = new XYLineChart_AWT_float(" Figure", "u(t)", "u(t)", u_t, 4);
		chart13.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart13 );          
		chart13.setVisible( true );

		XYLineChart_AWT_float chart14 = new XYLineChart_AWT_float(" Figure", "c(t)", "c(t)", c_t, 4);
		chart14.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart14 );          
		chart14.setVisible( true );

		XYLineChart_AWT_float chart15 = new XYLineChart_AWT_float(" Figure", "c_ut(t)", "c_ut(t)", c_ut_t, 4);
		chart15.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart15 );          
		chart15.setVisible( true );

		XYLineChart_AWT_float chart16 = new XYLineChart_AWT_float(" Figure", "i(t)", "i(t)", i_t, 4);
		chart16.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart16 );          
		chart16.setVisible( true );

		XYLineChart_AWT_float chart17 = new XYLineChart_AWT_float(" Figure", "a(t)", "a(t)", a_t, 4);
		chart17.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart17 );          
		chart17.setVisible( true );

		XYLineChart_AWT_float chart18 = new XYLineChart_AWT_float(" Figure", "h(t)", "h(t)", h_t, 4);
		chart18.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart18 );          
		chart18.setVisible( true );

		XYLineChart_AWT_float chart19 = new XYLineChart_AWT_float(" Figure", "rd(t)", "rd(t)", rd_t, 4);
		chart19.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart19 );          
		chart19.setVisible( true );

		XYLineChart_AWT_float chart20 = new XYLineChart_AWT_float(" Figure", "H(t)", "H(t)", H_t, 4);
		chart20.pack( );          
		RefineryUtilities.centerFrameOnScreen( chart20 );          
		chart20.setVisible( true );

		// XYLineChart_AWT_float chart21 = new XYLineChart_AWT_float(" Figure", "(t)", "i(t)", i_t, 4);
		// chart21.pack( );          
		// RefineryUtilities.centerFrameOnScreen( chart21 );          
		// chart21.setVisible( true );


	}

	public static class XYLineChart_AWT_float extends ApplicationFrame 
	{
	   public XYLineChart_AWT_float( String applicationTitle, String chartTitle, String name, float[] cseries, int len )
	   {
			super(applicationTitle);
			JFreeChart xylineChart = ChartFactory.createXYLineChart(
				chartTitle ,
				"Category" ,
				"Score" ,
				createDataset(name, cseries, len) ,
				PlotOrientation.VERTICAL ,
				true , true , false);
			 
			ChartPanel chartPanel = new ChartPanel( xylineChart );
			chartPanel.setPreferredSize( new java.awt.Dimension( 560 , 367 ) );
			final XYPlot plot = xylineChart.getXYPlot( );
			XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer( );
			renderer.setSeriesPaint( 0 , Color.RED );
			renderer.setSeriesPaint( 1 , Color.GREEN );
			renderer.setSeriesPaint( 2 , Color.YELLOW );
			renderer.setSeriesStroke( 0 , new BasicStroke( 4.0f ) );
			renderer.setSeriesStroke( 1 , new BasicStroke( 3.0f ) );
			renderer.setSeriesStroke( 2 , new BasicStroke( 2.0f ) );
			plot.setRenderer( renderer ); 
			setContentPane( chartPanel ); 
	   }
	   
	   public XYDataset createDataset( String name, float[] cseries, int len)
	   {
			final XYSeries series_Y = new XYSeries( name );
			for (int i = 0; i < len*T; i++){
				series_Y.add((double)(-80+i), (double)cseries[i]);
			}                   
			final XYSeriesCollection dataset = new XYSeriesCollection( );  
			dataset.addSeries( series_Y );
			return dataset;
	   }

	}

	public static class XYLineChart_AWT_long extends ApplicationFrame 
	{
	   public XYLineChart_AWT_long( String applicationTitle, String chartTitle, String name, long[] cseries, int len )
	   {
			super(applicationTitle);
			JFreeChart xylineChart = ChartFactory.createXYLineChart(
				chartTitle ,
				"Category" ,
				"Score" ,
				createDataset(name, cseries, len) ,
				PlotOrientation.VERTICAL ,
				true , true , false);
			 
			ChartPanel chartPanel = new ChartPanel( xylineChart );
			chartPanel.setPreferredSize( new java.awt.Dimension( 560 , 367 ) );
			final XYPlot plot = xylineChart.getXYPlot( );
			XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer( );
			renderer.setSeriesPaint( 0 , Color.RED );
			renderer.setSeriesPaint( 1 , Color.GREEN );
			renderer.setSeriesPaint( 2 , Color.YELLOW );
			renderer.setSeriesStroke( 0 , new BasicStroke( 4.0f ) );
			renderer.setSeriesStroke( 1 , new BasicStroke( 3.0f ) );
			renderer.setSeriesStroke( 2 , new BasicStroke( 2.0f ) );
			plot.setRenderer( renderer ); 
			setContentPane( chartPanel ); 
	   }
	   
	   public XYDataset createDataset( String name, long[] cseries, int len)
	   {
			final XYSeries series_Y = new XYSeries( name );
			for (int i = 0; i < len*T; i++){
				series_Y.add((double)(-80+i), (double)cseries[i]);
			}                   
			final XYSeriesCollection dataset = new XYSeriesCollection( );  
			dataset.addSeries( series_Y );
			return dataset;
	   }

	}

	public static class XYLineChart_AWT_int extends ApplicationFrame 
	{
	   public XYLineChart_AWT_int( String applicationTitle, String chartTitle, String name, int[] cseries, int len )
	   {
			super(applicationTitle);
			JFreeChart xylineChart = ChartFactory.createXYLineChart(
				chartTitle ,
				"Category" ,
				"Score" ,
				createDataset(name, cseries, len) ,
				PlotOrientation.VERTICAL ,
				true , true , false);
			 
			ChartPanel chartPanel = new ChartPanel( xylineChart );
			chartPanel.setPreferredSize( new java.awt.Dimension( 560 , 367 ) );
			final XYPlot plot = xylineChart.getXYPlot( );
			XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer( );
			renderer.setSeriesPaint( 0 , Color.RED );
			renderer.setSeriesPaint( 1 , Color.GREEN );
			renderer.setSeriesPaint( 2 , Color.YELLOW );
			renderer.setSeriesStroke( 0 , new BasicStroke( 4.0f ) );
			renderer.setSeriesStroke( 1 , new BasicStroke( 3.0f ) );
			renderer.setSeriesStroke( 2 , new BasicStroke( 2.0f ) );
			plot.setRenderer( renderer ); 
			setContentPane( chartPanel ); 
	   }
	   
	   public XYDataset createDataset( String name, int[] cseries, int len)
	   {
			final XYSeries series_Y = new XYSeries( name );
			for (int i = 0; i < len*T; i++){
				series_Y.add((double)(-80+i), (double)cseries[i]);
			}                   
			final XYSeriesCollection dataset = new XYSeriesCollection( );  
			dataset.addSeries( series_Y );
			return dataset;
	   }

	}
}