import tkinter as tk
from tkinter import Button, Canvas
from PIL import Image, ImageTk, ImageEnhance
import numpy as np
import matplotlib.pyplot as plt


class ImageCropperApp:
    def __init__(self, root, img1, img2, panels, brightness_factor, annotations_resized = None):
        self.root = root
        self.root.title("Image Cropper")

        self.panels = panels
        self.saved = False

        # Store original numpy images
        self.img1_8bit = img1
        self.img2_8bit = img2
        
        # Convert numpy arrays to PIL images
        self.img1 = self.array_to_image(img1)
        self.img2 = self.array_to_image(img2)

        ## Increase brightness manually before displaying 
        # To compensate the grayscale to RGB conversion for tkinter
        self.img1 = self.adjust_brightness(self.img1, brightness_factor)
        self.img2 = self.adjust_brightness(self.img2, brightness_factor)

        # Store original image sizes
        self.original_img1 = self.img1.copy()
        self.original_img2 = self.img2.copy()

        # Resize images for display
        self.canvas_width = 600
        self.canvas_height = 400
        self.img1_resized = self.resize_image(self.img1, self.canvas_width, self.canvas_height)
        self.img2_resized = self.resize_image(self.img2, self.canvas_width, self.canvas_height)

        # Create Canvases
        self.canvas1 = Canvas(root, width=self.canvas_width, height=self.canvas_height)
        self.canvas1.grid(row=0, column=0)
        self.canvas2 = Canvas(root, width=self.canvas_width, height=self.canvas_height)
        self.canvas2.grid(row=0, column=1)

        # Display images
        self.img1_tk = ImageTk.PhotoImage(self.img1_resized)
        self.img2_tk = ImageTk.PhotoImage(self.img2_resized)

        self.canvas1.create_image(0, 0, anchor=tk.NW, image=self.img1_tk)
        self.canvas2.create_image(0, 0, anchor=tk.NW, image=self.img2_tk)

        if annotations_resized:
            # Store annotations
            self.annotations_resized = annotations_resized

            # Draw annotations with scaling
            self.draw_annotations(self.canvas1, self.annotations_resized[0], "red", 
                                self.original_img1.width, self.original_img1.height, 
                                self.img1_resized.width, self.img1_resized.height)
            self.draw_annotations(self.canvas2, self.annotations_resized[1], "red", 
                                self.original_img2.width, self.original_img2.height, 
                                self.img2_resized.width, self.img2_resized.height)

        # Image Labels
        self.img1_label = tk.Label(root, text= f'Panel {self.panels[0]}')
        self.img1_label.grid(row=1, column=0)

        self.img2_label = tk.Label(root, text= f'Panel {self.panels[1]}')
        self.img2_label.grid(row=1, column=1)

        # Brightness sliders
        self.brightness_slider1 = tk.Scale(root, from_=0.1, to=40.0, resolution=0.1, orient=tk.HORIZONTAL,
                                        label="Brightness", command=self.update_brightness_img1)
        self.brightness_slider1.set(1.0)
        self.brightness_slider1.grid(row=5, column=0)

        self.brightness_slider2 = tk.Scale(root, from_=0.1, to=40.0, resolution=0.1, orient=tk.HORIZONTAL,
                                        label="Brightness", command=self.update_brightness_img2)
        self.brightness_slider2.set(1.0)
        self.brightness_slider2.grid(row=5, column=1)

        self.reset_button = Button(root, text="Reset Brightness", command=self.reset_brightness)
        self.reset_button.grid(row=5, column=0, columnspan=2)





        # Save status label
        self.save_status_label = tk.Label(root, text="", fg="green")
        self.save_status_label.grid(row=3, column=0, columnspan=2)

        # Variables for cropping
        self.crop_start1 = None
        self.crop_end1 = None
        self.rect1 = None
        self.crop_start2 = None
        self.crop_end2 = None
        self.rect2 = None

        # Buttons
        self.crop_button = Button(root, text="Crop", command=self.crop_image)
        self.crop_button.grid(row=2, column=0, columnspan=2)

        self.save_button = Button(root, text="Save Cropped Images", command=self.save_cropped_images)
        self.save_button.grid(row=4, column=0, columnspan=2)

        self.canvas1.bind("<ButtonPress-1>", self.on_press1)
        self.canvas1.bind("<B1-Motion>", self.on_drag1)
        self.canvas1.bind("<ButtonRelease-1>", self.on_release1)

        self.canvas2.bind("<ButtonPress-1>", self.on_press2)
        self.canvas2.bind("<B1-Motion>", self.on_drag2)
        self.canvas2.bind("<ButtonRelease-1>", self.on_release2)

        # Dictionary to store cropped images and their coordinates
        self.cropped_images = {"img1": None, "img2": None}

        # Add a Next button to close the window
        self.close_button = Button(root, text="Next", command=self.close_window)
        self.close_button.grid(row=6, column=0, columnspan=2)

    def array_to_image(self, img_array):
        """Convert numpy array to PIL image (assuming grayscale)."""
        return Image.fromarray(img_array).convert("L")
    
    def adjust_brightness(self, pil_image, factor):
        """Increase brightness of grayscale image while preserving mode."""
        enhancer = ImageEnhance.Brightness(pil_image)
        return enhancer.enhance(factor)  # Factor > 1 brightens, < 1 darkens

    def update_brightness_img1(self, event=None):
        """Update brightness of image 1 live based on slider."""
        self.brightness_factor1 = float(self.brightness_slider1.get())
        self.img1 = self.adjust_brightness(self.original_img1.copy(), self.brightness_factor1)
        self.img1_resized = self.resize_image(self.img1, self.canvas_width, self.canvas_height)
        self.img1_tk = ImageTk.PhotoImage(self.img1_resized)
        self.canvas1.delete("all")
        self.canvas1.create_image(0, 0, anchor=tk.NW, image=self.img1_tk)

        if hasattr(self, 'annotations_resized') and self.annotations_resized:
            # Draw annotations with scaling
            self.draw_annotations(self.canvas1, self.annotations_resized[0], "red", 
                                self.original_img1.width, self.original_img1.height, 
                                self.img1_resized.width, self.img1_resized.height)
            self.draw_annotations(self.canvas2, self.annotations_resized[1], "red", 
                                self.original_img2.width, self.original_img2.height, 
                                self.img2_resized.width, self.img2_resized.height)

    def update_brightness_img2(self, event=None):
        """Update brightness of image 2 live based on slider."""
        self.brightness_factor2 = float(self.brightness_slider2.get())
        self.img2 = self.adjust_brightness(self.original_img2.copy(), self.brightness_factor2)
        self.img2_resized = self.resize_image(self.img2, self.canvas_width, self.canvas_height)
        self.img2_tk = ImageTk.PhotoImage(self.img2_resized)
        self.canvas2.delete("all")
        self.canvas2.create_image(0, 0, anchor=tk.NW, image=self.img2_tk)

        if hasattr(self, 'annotations_resized') and self.annotations_resized:
            # Draw annotations with scaling
            self.draw_annotations(self.canvas1, self.annotations_resized[0], "red", 
                                self.original_img1.width, self.original_img1.height, 
                                self.img1_resized.width, self.img1_resized.height)
            self.draw_annotations(self.canvas2, self.annotations_resized[1], "red", 
                                self.original_img2.width, self.original_img2.height, 
                                self.img2_resized.width, self.img2_resized.height)
            
    def reset_brightness(self):
        self.brightness_slider1.set(1.0)
        self.brightness_slider2.set(1.0)
        self.update_brightness_img1(1.0)
        self.update_brightness_img2(1.0)

    def resize_image(self, img, width, height):
        img.thumbnail((width, height), Image.Resampling.LANCZOS)
        return img
    
    def scale_polygon(self, polygon, original_w, original_h, resized_w, resized_h):
        """Scale polygon coordinates to match resized image dimensions."""
        scale_x = resized_w / original_w
        scale_y = resized_h / original_h
        return [(x * scale_x, y * scale_y) for x, y in polygon.exterior.coords]

    def draw_annotations(self, canvas, polygons, color, original_w, original_h, resized_w, resized_h):
        """Draw scaled polygon annotations on the given canvas."""
        for polygon in polygons:
            scaled_coords = self.scale_polygon(polygon, original_w, original_h, resized_w, resized_h)
            coords_flattened = [coord for point in scaled_coords for coord in point]
            canvas.create_polygon(coords_flattened, outline=color, width=2, fill='', smooth=True)


    def on_press1(self, event):
        self.crop_start1 = (event.x, event.y)
        if self.rect1:
            self.canvas1.delete(self.rect1)
        self.rect1 = self.canvas1.create_rectangle(event.x, event.y, event.x, event.y, outline="red")

    def on_drag1(self, event):
        self.canvas1.coords(self.rect1, self.crop_start1[0], self.crop_start1[1], event.x, event.y)

    def on_release1(self, event):
        self.crop_end1 = (event.x, event.y)

    def on_press2(self, event):
        self.crop_start2 = (event.x, event.y)
        if self.rect2:
            self.canvas2.delete(self.rect2)
        self.rect2 = self.canvas2.create_rectangle(event.x, event.y, event.x, event.y, outline="red")

    def on_drag2(self, event):
        self.canvas2.coords(self.rect2, self.crop_start2[0], self.crop_start2[1], event.x, event.y)

    def on_release2(self, event):
        self.crop_end2 = (event.x, event.y)

    def crop_image(self):
        # Crop Image 1
        if self.crop_start1 and self.crop_end1:
            x1, y1 = self.crop_start1
            x2, y2 = self.crop_end1

            # Scale factor from displayed image to original
            scale_x = self.original_img1.width / self.img1_resized.width
            scale_y = self.original_img1.height / self.img1_resized.height

            # Convert to original coordinates
            orig_x1, orig_y1 = int(x1 * scale_x), int(y1 * scale_y)
            orig_x2, orig_y2 = int(x2 * scale_x), int(y2 * scale_y)

            # Crop the original np array
            cropped_img1 = self.img1_8bit[orig_y1:orig_y2, orig_x1:orig_x2]

            mean_val = cropped_img1.mean()
            brightened_cropped_img1 = cropped_img1 * self.brightness_factor1 + (1 - self.brightness_factor1) * mean_val
            # Clip values to stay within [0, 255]
            brightened_cropped_img1 = np.clip(brightened_cropped_img1, 0, 255).astype(np.float32)


            # Store both the cropped image and the crop coordinates
            self.cropped_images[self.panels[0]] = (brightened_cropped_img1, (orig_x1, orig_y1))

            # Update canvas
            self.img1 = self.array_to_image(cropped_img1)
            # Increase brightness manually before displaying to compensate the grayscale to RGB conversion for tkinter
            self.img1 = self.adjust_brightness(self.img1, self.brightness_factor1)

            self.img1_resized = self.resize_image(self.img1, self.canvas_width, self.canvas_height)
            self.img1_tk = ImageTk.PhotoImage(self.img1_resized)
            self.canvas1.delete("all")
            self.canvas1.create_image(0, 0, anchor=tk.NW, image=self.img1_tk)

        # Crop Image 2
        if self.crop_start2 and self.crop_end2:
            x1, y1 = self.crop_start2
            x2, y2 = self.crop_end2

            scale_x = self.original_img2.width / self.img2_resized.width
            scale_y = self.original_img2.height / self.img2_resized.height

            orig_x1, orig_y1 = int(x1 * scale_x), int(y1 * scale_y)
            orig_x2, orig_y2 = int(x2 * scale_x), int(y2 * scale_y)

            cropped_img2 = self.img2_8bit[orig_y1:orig_y2, orig_x1:orig_x2]

            mean_val = cropped_img2.mean()
            brightened_cropped_img2 = cropped_img2 * self.brightness_factor2 + (1 - self.brightness_factor2) * mean_val
            # Clip values to stay within [0, 255]
            brightened_cropped_img2 = np.clip(brightened_cropped_img2, 0, 255).astype(np.float32)

            self.cropped_images[self.panels[1]] = (brightened_cropped_img2, (orig_x1, orig_y1))

            self.img2 = self.array_to_image(cropped_img2)
            # Increase brightness manually before displaying to compensate the grayscale to RGB conversion for tkinter
            self.img2 = self.adjust_brightness(self.img2, self.brightness_factor2)

            self.img2_resized = self.resize_image(self.img2, self.canvas_width, self.canvas_height)
            self.img2_tk = ImageTk.PhotoImage(self.img2_resized)
            self.canvas2.delete("all")
            self.canvas2.create_image(0, 0, anchor=tk.NW, image=self.img2_tk)

    def save_cropped_images(self):
        # Here we simply update the label with the message
        self.save_status_label.config(text="Images Saved Successfully!")
        self.saved = True

    def close_window(self):

        cropped_img1 = self.img1_8bit
        mean_val = cropped_img1.mean()
        brightened_cropped_img1 = cropped_img1 * self.brightness_factor1 + (1 - self.brightness_factor1) * mean_val
        brightened_cropped_img1 = np.clip(brightened_cropped_img1, 0, 255).astype(np.float32)
        self.cropped_images[self.panels[0]] = (brightened_cropped_img1, (0, 0))

        cropped_img2 = self.img2_8bit
        mean_val = cropped_img2.mean()
        brightened_cropped_img2 = cropped_img2 * self.brightness_factor2 + (1 - self.brightness_factor2) * mean_val
        brightened_cropped_img2 = np.clip(brightened_cropped_img2, 0, 255).astype(np.float32)
        self.cropped_images[self.panels[1]] = (brightened_cropped_img2, (0, 0))

        self.root.quit() 
        self.root.destroy()
