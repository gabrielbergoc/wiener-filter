/*
 * To the extent possible under law, the ImageJ developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     http://creativecommons.org/publicdomain/zero/1.0/
 */

package com.mycompany.imagej;

import org.jtransforms.fft.FloatFFT_2D;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import lib.ImageAccess;

/**
 * A template for processing each pixel of either
 * GRAY8, GRAY16, GRAY32 or COLOR_RGB images.
 *
 * @author Johannes Schindelin
 */
public class Motion_Blur implements PlugInFilter {
	protected ImagePlus image;

	// plugin parameters
	private boolean direction;	// true: horizontal, false: vertical
	private String method;		// Convolution or Frequency multiplication
	private int filterSize;

	@Override
	public int setup(String arg, ImagePlus imp) {
		if (arg.equals("about")) {
			showAbout();
			return DONE;
		}

		image = imp;
		return DOES_8G | DOES_16 | DOES_32 | DOES_RGB;
	}

	@Override
	public void run(ImageProcessor ip) {
		if (showDialog()) {
			new ImageAccess(process(ip)).show("Result");
		}
	}

	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("Process pixels");

		gd.addChoice("Blur direction", new String[] { "Horizontal", "Vertical" }, "Horizontal");
		gd.addChoice("Method", new String[] { "Convolution", "Frequency multiplication" }, "Convolution");
		gd.addNumericField("Filter size", 9, 0);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		direction = gd.getNextChoice() == "Horizontal";
		method = gd.getNextChoice();
		filterSize = (int)gd.getNextNumber();

		return true;
	}

	public ImageProcessor process(ImageProcessor ip) {
		switch (method) {
			case "Convolution":
				return motionBlurConvolution(ip);

			case "Frequency multiplication":
				return motionBlurDFT(ip, filterSize);

			default:
				throw new Error("Invalid method");
		}
	}

	public ImageProcessor motionBlurConvolution(ImageProcessor ip) {
		ip.convolve(makeKernelConvolution(), direction ? 9 : 1, direction ? 1 : 9);
		return ip;
	}

	public ImageProcessor motionBlurDFT(ImageProcessor ip, int filterSize) {
		int nx = ip.getWidth();
		int ny = ip.getHeight();
		float[] kernel = makeKernelDFT(nx, ny, filterSize);
		float[] image = new float[2 * nx * ny];									// Array para conter a imagem no domínio da frequência
		float[] pixels = (float[]) ip.convertToFloatProcessor().getPixels();	// Array com a imagem no domínio espacial
		FloatFFT_2D fft = new FloatFFT_2D(ny, nx);								// Transformador FFT
		
		// Preencher array com informações no domínio espacial (depois da FFT,
		// ele conterá a imagem no domínio da frequência)
		for (int i = 0; i < nx * ny; i++) {
			image[2 * i] = pixels[i]; // Parte real
			image[2 * i + 1] = 0.0f; // Parte imaginária
		}

		// Executar a FFT na imagem e no kernel
		fft.complexForward(image);
		fft.complexForward(kernel);

		// Multiplicar as duas transformadas no domínio da frequência
		// (multiplicação de números complexos)
		float[] result = new float[2 * nx * ny];
		for (int i = 0; i < result.length; i += 2) {
			float re1 = image[i];
			float im1 = image[i + 1];
			float re2 = kernel[i];
			float im2 = kernel[i + 1];
			result[i] = re1 * re2 - im1 * im2;
			result[i + 1] = re1 * im2 + im1 * re2;
		}

		// Executar a transformada inversa para voltar ao domínio espacial
		fft.complexInverse(result, true);

		// Atualizar o array do domínio espacial
		for (int i = 0; i < nx * ny; i++) {
			pixels[i] = result[2 * i];
		}

		return new FloatProcessor(nx, ny, pixels);
	}

	private float[] makeKernelConvolution() {
		return makeKernelConvolution(9);
	}

	private float[] makeKernelConvolution(int size) {
		float[] kernel = new float[size];

		for (int i = 0; i < size; i++) {
			kernel[i] = 1f / size;
		}

		return kernel;
	}

	private float[] makeKernelDFT(int width, int height, int size) {
		float[] kernel = new float[2 * width * height]; // Espaço para parte real e imaginária

		// Popular parte real do kernel (parte imaginária é sempre zero)
		for (int i = 0; i < width * height; i++) {
			if (direction) {
				kernel[2 * i] = i < size ? 1.0f / size : 0;
			} else {
				kernel[2 * i] = i % width == 0 && i / width < size ? 1.0f / size : 0;
			}
		}

		return kernel;
	}

	public void showAbout() {
		IJ.showMessage("MotionBlur",
			"a motion blur filter - works horizontally and vertically, by convolution or frequency multiplication");
	}

	/**
	 * Main method for debugging.
	 *
	 * For debugging, it is convenient to have a method that starts ImageJ, loads
	 * an image and calls the plugin, e.g. after setting breakpoints.
	 *
	 * @param args unused
	 */
	public static void main(String[] args) throws Exception {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		// see: https://stackoverflow.com/a/7060464/1207769
		Class<?> clazz = Motion_Blur.class;
		java.net.URL url = clazz.getProtectionDomain().getCodeSource().getLocation();
		java.io.File file = new java.io.File(url.toURI());
		System.setProperty("plugins.dir", file.getAbsolutePath());

		// start ImageJ
		new ImageJ();

		// open the Clown sample
		ImagePlus image = IJ.openImage("http://imagej.net/images/clown.jpg");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
