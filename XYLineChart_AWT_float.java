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


public class XYLineChart_AWT_float extends ApplicationFrame
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
		chartPanel.setPreferredSize( new java.awt.Dimension( 660 , 467 ) );
		final XYPlot plot = xylineChart.getXYPlot( );
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer( );
		renderer.setSeriesPaint( 0 , Color.BLACK );
		renderer.setSeriesPaint( 1 , Color.RED );
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
		for (int i = 0; i < len*TheoryModel.T; i++){
			series_Y.add((double)(-80+i), (double)cseries[i]);
		}                   
		final XYSeriesCollection dataset = new XYSeriesCollection( );  
		dataset.addSeries( series_Y );
		return dataset;
   }

}