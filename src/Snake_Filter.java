/**
 *
 * @author  Javier
 * @author  Nicolas
 * @author  Andy
 * @author  Harvy
 * @author  Mileidy
 * @author  Maria Alejandra
 * @author  Juan David
 */
import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.Toolbar;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class Snake_Filter implements PlugIn {

    /**
     *imp
     */
    public ImagePlus imp; 	
        
	/** Procesador actualmente activo */
	public ImageProcessor ip;
	
        /** Radio de suavizado utilizado para calcular el gradiente **/
	float dSuavizado;
	
        /** Derivadas de gaussiano (primer orden en x e y) **/
	ImageProcessor gaussianX;
	ImageProcessor gaussianY;
        
        /** Campo vectorial normalizado de gradiente **/
	float [][][] campoNormalizado;
	
        /** Paso de tiempo de integración **/
	float dTiempoEjecucion;
	
	/** Cuadrícula de espacios ocupados **/
	List[][] Grid; // = new List<String>[2]; Inicializar la variable
	
	/** Tamaño de la cuadrícula **/
	int nWG,nHG;

	/** nuevas coordenadas de puntos de semilla **/
	ArrayList<Float []> puntosBase = new ArrayList<>(); 
	
	/** distancia de separación entre curvas de nivel **/
	float distanciaSeparacion;
        
	/** Detener la integración de la línea cuando está más cerca que * (dPararLineaSeparacion * distanciaSeparacion) **/
	float dPararLineaSeparacion;
	
	/** Ignorar los contornos abiertos **/
	boolean bEliminarAbierto;
        
        /** Tamaños de la grilla **/
	float dGridSize;
	
        /** Tamaño de la imagen **/
	int nW,nH;

	/** nombre de LUT elegido **/
	String sNombreLut;

	/** Matriz de colores LUT **/	
	int [][] RGBLutTabla;
	
	/** invertir LUT **/
	boolean bInvertirLUT;
	
	/** si usar un solo color o LUT **/
	boolean bColorSimple;
	
	/** color seleccionado de Selector de color / Barra de herramientas **/
	Color sistemaColor;
	
	//double dStrokeWidth;  // Anchura del trazo
	
	/** superposición para renderizar **/
	Overlay image_overlay; 

	/** Intensidades mínimas y máximas de la imagen tomadas de Contraste / Brillo **/
	float fInMin;
	float fInRango;
	float fInMax;
	
	/** semilla inicial para la primera línea,
	 * punto con la pendiente mínima **/
	Point iniSemilla;
	
	/** número total de contornos **/
	int nLineas;
	
	@Override
	public void run(String arg) {
		
		PolygonRoi new_line;
		
		/** variables de medición del tiempo **/
		long startTime, contourTime=0;
		
		imp = IJ.getImage();
		if(null == imp)
		{
                    IJ.noImage();
                    return;
		}
		else if (imp.getType() != ImagePlus.GRAY8 && imp.getType() != ImagePlus.GRAY16 && imp.getType() != ImagePlus.GRAY32) 
		{
		    IJ.error("8, 16 or 32-bit imagen en escala de grises requerida");
		    return;
		}	
		// pedir parámetros

		if(!showParametersDialog())
                    return;
		
		IJ.log("Filtro Snake ( contornos )");
		IJ.log("Parametros:");
		IJ.log("Radio de suavizado: "+Float.toString(dSuavizado)+" pixels");
		IJ.log("Paso de integración de linea: "+Float.toString(dTiempoEjecucion)+" pixels");
		IJ.log("Distancia entre lineas: "+Float.toString(distanciaSeparacion)+" pixels");
		IJ.log("Linea de parada a una fracción de la distancia: "+Float.toString(dPararLineaSeparacion));
		
		if(bColorSimple)
                    IJ.log("Usando un solo color");
		else
		{
                    if(bInvertirLUT)
                        IJ.log("Utilizando "+ sNombreLut +" LUT (Invertido)");
                    else
                        IJ.log("Utilizando "+ sNombreLut);
		}
		if(bEliminarAbierto)
		{
			IJ.log("Se eliminan los contornos abiertos");
		}
		// Comencemos a medir el tiempo
		startTime = System.nanoTime();
				
		nW=imp.getWidth();
		nH=imp.getHeight();
		ip=imp.getProcessor();

		/** Tamaño de la cuadrícula **/
		dGridSize=0.5f*distanciaSeparacion;
		nWG=(int)Math.ceil(nW/(float)dGridSize);
		nHG=(int)Math.ceil(nH/(float)dGridSize);

		Grid = new ArrayList[nWG][nHG];

		/** Derivadas del cálculo gaussiano **/
		double[][] kernel = computeKernelGaussian2D_du(dSuavizado,dSuavizado, 0.0f);		
		gaussianX= convolve(ip, kernel);
		kernel = computeKernelGaussian2D_dv(dSuavizado,dSuavizado, 0.0f);
                gaussianY= convolve(ip, kernel);		
		
		// normalizar el campo de gradiente
		campoNormalizado  = new float [nW][nH][2];
		fillNormalizedField();
		image_overlay = new Overlay();
		
		// configurar los colores de los contornos

		if (bColorSimple)
                    sistemaColor = Toolbar.getForegroundColor();
		else
		{
                    getRGBLutTable();
                    fInMin = (float) ip.getMin();
                    fInMax=(float) ip.getMax();
                    fInRango = fInMax -fInMin;
		}
		
		nLineas=0;
                
		// Primera linea
		new_line=buildLine(iniSemilla.getX(),iniSemilla.getY());
		if(new_line!=null)
		{
                // Agregar a la superposición
                image_overlay.add(new_line);		
                nLineas++;
                IJ.log("Linea "+ Integer.toString(nLineas)+" Añadir.");
		}

		// Abrir la ventana
		imp.setOverlay(image_overlay);
		imp.updateAndRepaintWindow();
		imp.show();
		
		// Utilizando semillas de la primera línea de contorno

		if (generateLinesFromSeedList())
		{
                    // Comprobar semillas de celdas vacías
                    refillSeedList();
                    // Reintento
                    generateLinesFromSeedList();
		}
		
		IJ.log("Hecho");
		contourTime = System.nanoTime() - startTime;
		IJ.log("Tiempo Total: " + String.format("%.2f",((double)Math.abs(contourTime))*0.000000001) + " s");
	}
	
	/** Generar líneas usando la cola puntosBase **/
	boolean generateLinesFromSeedList()
	{
            // Pasando por la cola de puntos semilla y Generando por lineas
            PolygonRoi new_line;
            boolean goOn=true;
            boolean bFoundSeed; 
            Float [] seed_point = new Float [2];
            while (goOn)
            {
                if(puntosBase.isEmpty())
                {
                    goOn=false;
                }
                else
                {
                    bFoundSeed = false;

                    while (!bFoundSeed)
                    {
                        // La cola está vacía, terminar
                        if(puntosBase.isEmpty())
                        {
                            goOn = false;
                            break;
                        }
                        seed_point= puntosBase.get(0);

                        if(isBusy(seed_point[0],seed_point[1],0.9f*distanciaSeparacion))
                        {
                            // Si el punto de semilla está demasiado cerca de otra línea de contorno,saca la semilla de la cola
                            puntosBase.remove(0);
                        }
                        else
                        {
                            // El punto de semilla está bien
                            bFoundSeed=true;				
                            puntosBase.remove(0);
                        }
                    }
                    
                    // Generar una línea de contorno desde el punto inicial
                    if(goOn && bFoundSeed)
                    {
                        new_line=buildLine(seed_point[0],seed_point[1]);
                        if (new_line!=null)
                        {
                            if(new_line.size()>2)
                            {
                                image_overlay.add(new_line);
                                nLineas++;
                                imp.setOverlay(image_overlay);
                                imp.updateAndRepaintWindow();
                                imp.show();
                                IJ.log("Linea "+ Integer.toString(nLineas)+" Añadir.");
                                
                                if(!imp.isProcessor())
                                {
                                    IJ.log("Imagen cerrada. Terminando.");
                                    return false;
                                }
                            }
                        }
                    }
                }
            }
            return true;

	}
	/** Rutina de integración de línea utiliza el método de Newton simple con un paso definido por el usuario (paso de integración)
            * @param xstart
            * @param ystart
            * @return  **/
	public PolygonRoi buildLine(float xstart, float ystart)
	{
		/** Propia cuadrícula de la línea actual para verificar si hay auto-intersecciones.
		  Es igual al tamaño de la imagen. La precisión de más de 1 píxel parece innecesaria **/
		int [][] ownGrid  = new int [nW][nH];
		int gx, gy, gxn,gyn,nDirection;
                
                PolygonRoi polyline;
                
		float xcurr, ycurr, xvel, yvel;
		ArrayList<Float> px = new ArrayList<>(); 
		ArrayList<Float> py = new ArrayList<>();

		boolean bEnd;
		
		px.add(xstart);
		py.add(ystart);
		
		// La línea crecerá en dirección hacia adelante y hacia atrás.
		//  Esta variable lo especifica.
		nDirection =1;
		
		for(nDirection=1;nDirection>-2;nDirection-=2)
		{
			bEnd = false;
			xcurr=xstart;
			ycurr=ystart;
                        
			// Marcar el píxel OwnGrid como actual (1)
			gx = (int) Math.floor(xstart);
			gy = (int) Math.floor(ystart);
			ownGrid[gx][gy]=1;

			while (!bEnd)
			{
                            // Gradiente en la posición actual (+/- depende de nDirection)
                            xvel = ((float)nDirection)*campoNormalizado[Math.round(xcurr)][Math.round(ycurr)][0];
                            yvel = ((float)nDirection)*campoNormalizado[Math.round(xcurr)][Math.round(ycurr)][1];

                            // Compruebe que no sea un sumidero o una fuente o simplemente cero píxeles
                            if(Float.isNaN((Float)xvel))
                            {
                                bEnd=true;
                                ownGrid[gx][gy]=2;
                            }
                            else
                            {
                                if(Float.isNaN((Float)xcurr))
                                {
                                    bEnd = true;
                                    ownGrid[gx][gy]=2;
                                }
                                else
                                {
                                    // Hacer crecer una línea por paso de integración
                                    xcurr = xcurr +xvel * dTiempoEjecucion;
                                    ycurr = ycurr +yvel * dTiempoEjecucion;

                                    // La línea está cerca de otra línea, abortar el crecimiento
                                    if(isBusy(xcurr,ycurr,distanciaSeparacion*dPararLineaSeparacion))
                                    {
                                        bEnd = true;
                                        ownGrid[gx][gy]=2;
                                    }
                                    else
                                    {
                                        // Revisemos la auto-intersección
                                        gxn = (int) Math.floor(xcurr);
                                        gyn = (int) Math.floor(ycurr);
                                        // Fuera de imagen, abortar!
                                        if(gxn<0 || gyn<0 || gxn>(nW-1)|| gyn>(nH-1))
                                        {
                                            bEnd = true;
                                        }
                                        else
                                        {
                                            // auto-intersección, abortar!
                                            if(ownGrid[gxn][gyn]==2)
                                            {
                                                bEnd = true;
                                            }
                                            else
                                            {
                                                // Salimos del píxel de "crecimiento" actual
                                                // Vamos a marcarlo con (2) = visitado
                                                // y marcar el nuevo píxel como actual (1)
                                                if(gxn!=gx || gyn!=gy)
                                                {
                                                    ownGrid[gx][gy]=2;
                                                    gx=gxn;
                                                    gy=gyn;
                                                    ownGrid[gx][gy]=1;
                                                }
                                                // agregar un punto al final de la matriz de líneas (crecimiento hacia adelante)
                                                if(nDirection==1)
                                                {
                                                    px.add(xcurr);
                                                    py.add(ycurr);

                                                }
                                                // Agregar un punto al comienzo de la matriz de líneas (crecimiento hacia atrás)
                                                // de esta manera los puntos en la matriz se indexan continuamente
                                                else
                                                {
                                                        px.add(0,xcurr);
                                                        py.add(0,ycurr);								
                                                }
                                            }
                                        }
                                    }
                                }
                            }   
                        }
                // Estado aquí también
                ownGrid[gx][gy]=2;
            }
		//Variables de generacion linea
		float[] floatX = new float[px.size()];
		float[] floatY = new float[px.size()];
		float nAverVal=0;
		int i = 0;
		int nPoints =px.size();
		for (i=0;i<nPoints;i++) 
		{
                    // Conviértalo en dos matrices que podrían alimentarse a PolygonRoi
                    floatX[i]=px.get(i);
                    floatY[i]=py.get(i);
                    
                    // Agregue puntos semilla en ambos lados de cada punto en la línea
                    addSeedPoint(floatX[i],floatY[i]);
                    // Agregue una nueva línea a la cuadrícula que rastrea el espacio ocupado
                    markOccupied(px.get(i), py.get(i));
		}
		
		// ROI de la línea de contorno		
		// Vamos a comprobar la distancia desde el principio hasta el final
		// Si tiene menos de 1,5 píxeles, conviértalo en un contorno cerrado

		 if(Math.sqrt(Math.pow(floatX[0]-floatX[nPoints-1], 2)+Math.pow(floatY[0]-floatY[nPoints-1], 2))<1.5)
		 {
                    polyline = new PolygonRoi(floatX, floatY, Roi.POLYGON);
		 }
		 else //polyline
		 {
                    if(bEliminarAbierto)
                    {	 return null;}
                    else
                    {polyline = new PolygonRoi(floatX, floatY, Roi.POLYLINE);}
		 }
                 
		 // Un solo color
		 if(bColorSimple)
		 {
                    polyline.setStrokeColor(sistemaColor);
		 }
                 
		 // Uso de colores LUT para codificar la profundidad del color
		 else
		 {
                    // Calcular la intensidad media por línea
                    nAverVal=0;
                    
                    for (i=0;i<nPoints;i++) 
                    {
                            nAverVal += ip.getf((int)Math.floor(floatX[i]), (int)Math.floor(floatY[i]));
                    }
                    nAverVal/=px.size();	

                    // Elegir la configuración actual de contraste / brillo de la imagen
                    if(nAverVal>fInMax)
                            nAverVal=255;
                    else
                    {
                            if(nAverVal<fInMin)
                                nAverVal=0;
                            else
                                nAverVal=Math.round(255.0f*(nAverVal-fInMin)/fInRango);
                    }

                    polyline.setStrokeColor(new Color(RGBLutTabla[(int)nAverVal][0],RGBLutTabla[(int)nAverVal][1],RGBLutTabla[(int)nAverVal][2]));
		 }
		 
		//polyline.setStrokeWidth(dStrokeWidth);
		
		return polyline;
		
		
	}
	/** La función agrega un punto a la variable Grid
            * que rastrean la densidad local de las curvas de nivel en la imagen **/
	void markOccupied(float xin, float yin)
	{
		int gx,gy;
		
		ArrayList<Point> currarr;		
		
		//grid's cell coordinate
		gx = (int) Math.floor(xin/dGridSize);
		gy = (int) Math.floor(yin/dGridSize);
		currarr =(ArrayList<Point>) Grid[gx][gy];
		if(currarr==null)
		{
			Grid[gx][gy] = new ArrayList<Point>();
			currarr =(ArrayList<Point>) Grid[gx][gy];
		}
		currarr.add(new Point(xin,yin));
	}
	
	/** La función agrega dos puntos de semilla a la cola puntosBase
            (si no están en una forma de otra línea, es decir, más lejos
            que el valor distanciaSeparacion) **/
	void addSeedPoint(float xs, float ys)
	{
            float xn,yn,xvel, yvel;
            Float[] spoint;
            int nDirection;
            nDirection =1;

            for(nDirection=1;nDirection>-2;nDirection-=2)
            {
                // Dirección perpendicular
                xvel = campoNormalizado[Math.round(xs)][Math.round(ys)][1]*((float)nDirection);
                yvel = -campoNormalizado[Math.round(xs)][Math.round(ys)][0]*((float)nDirection);

                // Punto especial (plano o fuente o fregadero)
                if(Float.isNaN((Float)xvel) ||Float.isNaN((Float)yvel))
                {
                    return;
                }
                xn=xs+xvel*distanciaSeparacion*1.1f;
                yn=ys+yvel*distanciaSeparacion*1.1f;

                if(!(xn<0 || yn<0 || xn>(nW-1)|| yn>(nH-1)))
                {
                    spoint = new Float [2];
                    spoint[0]=xn;
                    spoint[1]=yn;
                    puntosBase.add(spoint);
                }
            }		
	}
	
	/** La función genera puntos de semilla adicionales en las celdas de Grid.
            * donde no hay puntos (para llenar la imagen completa **/
	void refillSeedList()
	{
            int gx,gy;
            ArrayList<Point> currarr;		
            Float[] spoint;

            //grid's cell coordinate
            for (gx=0;gx<nWG;gx++)
                for (gy=0;gy<nHG;gy++)
                {
                    currarr =(ArrayList<Point>) Grid[gx][gy];
                    //Limpiar celdas
                    if(currarr==null)
                    {
                        // Agregue un punto de semilla en el medio de la celda
                        spoint = new Float [2];
                        spoint[0]=((float)gx)*dGridSize+0.5f*dGridSize;
                        spoint[1]=((float)gy)*dGridSize+0.5f*dGridSize;
                        puntosBase.add(spoint);
                    }
                }
            }
	/** función comprobando si el punto xc, yc está cerca de cualquier otro existente
           * puntos (almacenados en Grid) por la distancia cDist
           * @param xc
           * @param yc
           * @param dDist
           * @return  **/
	public boolean isBusy(float xc, float yc, float dDist)
	{
		int gx, gy, x, y, i, j, k;
		ArrayList<Point> currarr;
		Point point;
		int nRange=(int) Math.ceil(dDist/dGridSize);
		
		x= Math.round(xc);
		y= Math.round(yc);
                
		// fuera de imagen
		if(x<0 || y<0 || x>(nW-1)|| y>(nH-1))
		{
                    return true;
		}
                
		// índice de celda de cuadrícula
		gx = (int) Math.floor(xc/dGridSize);
		gy = (int) Math.floor(yc/dGridSize);
		
		// Comprobando las celdas vecinas también
		for (i=gx-nRange;i<=gx+nRange;i++)
                    for (j=gy-nRange;j<=gy+nRange;j++)
                    {
                        if(!(i<0 || j<0 || i>(nWG-1)|| j>(nHG-1)) )
                        {
                            currarr =(ArrayList<Point>) Grid[i][j];
                            if(!(currarr==null))
                            {
                                for (k=0;k<currarr.size();k++)
                                    {
                                        point= currarr.get(k);
                                        if(point.distance(xc, yc)<dDist)
                                                return true;
                                    }
                            }
                        }
                    }
	
		return false;
		
	}
	/**  Calcula el campo vectorial normalizado que es ortogonal al campo de gradiente 
            * calculado en los procesadores de imágenes gaussianX y gaussianY.
            * Se almacena en la variable campoNormalizado.
            * Además, almacena el punto con el gradiente de valor absoluto más pequeño en la variable iniSemilla **/
	public void fillNormalizedField()
	{
		int i,j,imin,jmin;
		float dx, dy, len, minlen=Float.MAX_VALUE;
		imin=-1;
		jmin=-1;
		for (i=0;i<nW;i++)
                    for (j=0;j<nH;j++)
                    {
                        dx = gaussianX.getf(i, j);
                        dy = gaussianY.getf(i, j);
                        len = (float) Math.sqrt(dx*dx+dy*dy);
                        if(len>0.00000001f)
                        {
                            campoNormalizado[i][j][0]=dy/len;
                            campoNormalizado[i][j][1]=-dx/len;
                            if (len<minlen)
                            {
                                minlen=len;
                                imin=i;
                                jmin=j;
                            }
                        }
                        else
                        {
                            campoNormalizado[i][j][0]=Float.NaN;
                            campoNormalizado[i][j][1]=Float.NaN;					
                        }
                    }
		// caso trivial, la imagen es completamente incorrecta
		// escojamos el punto medio de la imagen
		if(imin==-1)
		{
                    iniSemilla=new Point((float)nW*0.5f,(float)nH*0.5f);
		}
		else
		{
                    iniSemilla=new Point((float)imin,(float)jmin);
		}
		
	}
	/** Diálogo con parámetros de vinculación
            * @return  **/
	public boolean showParametersDialog()
	{
		int nLutChoice;
		GenericDialog contourlinesD = new GenericDialog("Rendering parameters");
		String [] luts = IJ.getLuts();
		
		contourlinesD.addNumericField("Radio de suavizado:", Prefs.get("ContourLines.dSmoothR", 2.0), 1, 3,"pixels");
		contourlinesD.addNumericField("Paso de integracion de linea (0.01-1):", Prefs.get("ContourLines.dTimeStep", 0.5), 2, 4,"pixels");
		contourlinesD.addNumericField("Distancia entre lineas:", Prefs.get("ContourLines.dSep", 5), 1, 3,"pixels");
		contourlinesD.addNumericField("Linea de parada a una fraccion de la distancia:", Prefs.get("ContourLines.dStopLineSep", 0.5), 1, 3,"fraction");
		contourlinesD.addCheckbox("Usar un solo color (actual) para dibujar lineas?", Prefs.get("ContourLines.bSingleColor", true));
		contourlinesD.addChoice("Contornos de codigo de color con LUT:",luts,Prefs.get("ContourLines.sLutChoice","Fire"));
		contourlinesD.addCheckbox("Invertir LUT?", Prefs.get("ContourLines.bInvertLUT", false));
		contourlinesD.addCheckbox("Eliminar contornos abiertos?", Prefs.get("ContourLines.bRemoveOpen", false));
		contourlinesD.setResizable(false);
		contourlinesD.showDialog();	
		if (contourlinesD.wasCanceled())
            return false;

		dSuavizado = (float) contourlinesD.getNextNumber();
		Prefs.set("ContourLines.dSmoothR", dSuavizado);
		dTiempoEjecucion = (float) contourlinesD.getNextNumber();
		Prefs.set("ContourLines.dTimeStep", dTiempoEjecucion);
		distanciaSeparacion = (float) contourlinesD.getNextNumber();
		Prefs.set("ContourLines.dSep", distanciaSeparacion);
		dPararLineaSeparacion = (float) contourlinesD.getNextNumber();
		Prefs.set("ContourLines.dStopLineSep", dPararLineaSeparacion);

		bColorSimple = contourlinesD.getNextBoolean();
		Prefs.set("ContourLines.bSingleColor", bColorSimple);
		nLutChoice = contourlinesD.getNextChoiceIndex();
		Prefs.set("ContourLines.sLutChoice", luts[nLutChoice]);
		sNombreLut = luts[nLutChoice];
		bInvertirLUT = contourlinesD.getNextBoolean();
		Prefs.set("ContourLines.bInvertLUT", bInvertirLUT);
		bEliminarAbierto = contourlinesD.getNextBoolean();
		Prefs.set("ContourLines.bRemoveOpen", bEliminarAbierto);

		return true;
	}
	/** Derivado 2D del kernel normalizado gaussiano en x (tomado del código Jalmar, podría ser excesivo)
            * @param sigma_x
            * @param sigma_y
            * @param theta
            * @return  **/
	public static double[][] computeKernelGaussian2D_du(double sigma_x, double sigma_y, float theta)
	{
            // calcular el tamaño de kernel requerido (2 * 3 * sigma ~ 99%)
            int kernel_radius = (int)Math.round(3*Math.max(sigma_x, sigma_y)); // RSLV: use floor instead of round?
            int kernel_size = 1+2*kernel_radius;

            // computar kernel
            double[][] kernel = new double[kernel_size][kernel_size];
            for(int ky = 0; ky < kernel_size; ++ky)
            {
                int y = ky - kernel_radius;
                for(int kx = 0; kx < kernel_size; ++kx)
                {
                    int x = kx - kernel_radius;
                    double u = x * Math.cos(theta) - y * Math.sin(theta);
                    double v = x * Math.sin(theta) + y * Math.cos(theta);
                    kernel[kx][ky] = gaussian2D_dx(u, v, sigma_x, sigma_y);
                }
            }
		
            // normalizar el kernel
            kernel = normalize_kernel(kernel);

            return kernel;
	}
	/** Derivado 2D del kernel normalizado gaussiano en y (tomado del código Jalmar, podría ser excesivo)
            * @param sigma_x
            * @param sigma_y
            * @param theta
            * @return  **/
	public static double[][] computeKernelGaussian2D_dv(double sigma_x, double sigma_y, float theta)
	{
            // calcular el tamaño de kernel requerido (2 * 3 * sigma ~ 99%)
            int kernel_radius = (int)Math.round(3*Math.max(sigma_x, sigma_y)); // RSLV: ¿usar piso en lugar de redonda?
            int kernel_size = 1+2*kernel_radius;

            // Computar kernel
            double[][] kernel = new double[kernel_size][kernel_size];
            for(int ky = 0; ky < kernel_size; ++ky)
            {
                int y = ky - kernel_radius;
                for(int kx = 0; kx < kernel_size; ++kx)
                {
                    int x = kx - kernel_radius;
                    double u = x * Math.cos(theta) - y * Math.sin(theta);
                    double v = x * Math.sin(theta) + y * Math.cos(theta);
                    kernel[kx][ky] = gaussian2D_dy(u, v, sigma_x, sigma_y);
                }
            }

            // normalizar el kernel
            kernel = normalize_kernel(kernel);

            return kernel;
	}
	/** Derivada de la función gaussiana en x
            * @param x
            * @param y
            * @param sigma_x
            * @param sigma_y
            * @return  **/
	public static double gaussian2D_dx(double x, double y, double sigma_x, double sigma_y)
	{
            return ((-x)/(2*Math.PI*Math.pow(sigma_x, 3)*sigma_y))*Math.exp(-0.5*((x*x)/(sigma_x*sigma_x)+(y*y)/(sigma_y*sigma_y)));
	}
	/** Derivada de la función gaussiana en y
            * @param x
            * @param y
            * @param sigma_x
            * @param sigma_y
            * @return  **/
	public static double gaussian2D_dy(double x, double y, double sigma_x, double sigma_y)
	{
            return ((-y)/(2*Math.PI*sigma_x*Math.pow(sigma_y, 3)))*Math.exp(-0.5*((x*x)/(sigma_x*sigma_x)+(y*y)/(sigma_y*sigma_y)));
	}
	/** Normalización del kernel DoG
        * @param kernel
        * @return  **/
	public static double[][] normalize_kernel(double[][] kernel)
	{
            // Calcular la suma de componentes
            double sum = 0.0;
            for(int kx = 0; kx < kernel.length; ++kx)
            {
                for(int ky = 0; ky < kernel[kx].length; ++ky)
                {
                    sum += Math.abs(kernel[kx][ky]); // NOTA: use abs para normalizar el núcleo simétrico con un lóbulo positivo y negativo
                }
            }
            // Evitar la división por cero
            if(sum == 0.0) { return kernel; }

            // Calcular el factor de escala
            double scale_factor = 1 / sum;

            // Componentes de escala
            for(int kx = 0; kx < kernel.length; ++kx)
            {
                for(int ky = 0; ky < kernel[kx].length; ++ky)
                {
                    kernel[kx][ky] *= scale_factor;
                }
            }

            return kernel;
	}
	/** Convolución de imagen con kernel DoG
        * @param ip
        * @param kernel
        * @return  **/
	public static ImageProcessor convolve(ImageProcessor ip, double[][] kernel)
	{
            // Obtener tamaños de imagen y kernel
            int image_width = ip.getWidth(),
                image_height = ip.getHeight(), 
                kernel_width = kernel.length,  
                kernel_height = kernel_width, // NOTA: suponga núcleo cuadrado
                kernel_half_width;
                kernel_half_width = (int)Math.floor(0.5 * kernel_width);
            int kernel_half_height = kernel_half_width;

            // Convertir el procesador de imagen de entrada a flotar
            ImageProcessor ip_inp = ip.convertToFloat(); // RSLV: duplicate first?

            // Crear un nuevo procesador flotante de salida vacío
            ImageProcessor ip_res = new FloatProcessor(image_width, image_height);

            // Convolucionar imagen de entrada con kernel
            for(int py = 0; py < image_height; ++py)
            {
                for(int px = 0; px < image_width; ++px)
                {
                    double kernel_product = 0.0;
                    for(int ky = 0; ky < kernel_height; ++ky)
                    {
                        int ppy = py + ky - kernel_half_height;
                        if(ppy < 0) ppy = 0; // Abrazadera en el borde 
                        if(ppy >= image_height) ppy = image_height - 1; // Abrazadera en el borde
                        for(int kx = 0; kx < kernel_width; ++kx)
                        {
                            int ppx = px + kx - kernel_half_width;
                            if(ppx < 0) ppx = 0; // Abrazadera en el borde 
                            if(ppx >= image_width) ppx = image_width - 1; // Abrazadera en el borde
                            kernel_product += ip_inp.getf(ppx, ppy) * kernel[kx][ky];
                        }
                    }
                    ip_res.setf(px, py, (float)kernel_product);
                }
            }
            return ip_res;
	}
	
	/** la función obtiene la LUT especificada por sZLUTName en la configuración
	 * y devuelve un mapa de tabla de 256x3 en formato RGB.
	 * Esto se hace principalmente desde que el LutLoader de ImageJ genera
	 * algunos LUT sobre la marcha (fuego, etc.) y algunos están cargados
	 * desde el disco duro. Además, aplica automáticamente LUT a la imagen actual.
	 * Entonces, para obtener una tabla de LUT, hago una imagen ficticia,
	 * codifíquelo con todos los colores de LUT y léalo en color de píxeles (y cierre la imagen).
	 * Este truco está tomado del complemento DoM */
	void  getRGBLutTable()
	{
		int i, j;
		int [] onepix; 
		RGBLutTabla = new int[256][3];
		ByteProcessor ish = new ByteProcessor(256,10);
		for (i=0; i<256; i++)
                    for (j=0; j<10; j++)
                        ish.putPixel(i, j, i);
		ImagePlus imLUT = new ImagePlus("Prueba",ish);
		imLUT.show();
		IJ.run(sNombreLut);
		IJ.run("RGB Color");

		imLUT.setSlice(1);
		for(i=0;i<256;i++)
		{			
                    onepix= imLUT.getPixel(i, 2);		
                    if(!bInvertirLUT)
                    {
                        for (j=0;j<3;j++)
                            {RGBLutTabla[i][j]=onepix[j];}
                    }
                    else
                    {
                        for (j=0;j<3;j++)
                            {RGBLutTabla[255-i][j]=onepix[j];}
                    }
		}
            imLUT.changes=false;
            imLUT.close();
	}
}