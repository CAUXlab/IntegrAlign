import numpy as np
import tkinter as tk
from tkinter import Button
from PIL import Image, ImageTk, ImageEnhance
from scipy.ndimage import rotate, shift
import gc
import matplotlib.pyplot as plt

class ImageManualAlignmentApp:
    def __init__(self, master, img1_8bit, img2_8bit, panels, brightness_factor):
        self.master = master
        self.master.title("Image Alignment Tool")

        self.panels = panels

        # Store original images as grayscale
        self.img1_original = Image.fromarray(img1_8bit).convert("L")  # Keep as grayscale
        self.img2_original = Image.fromarray(img2_8bit).convert("L")  # Keep as grayscale

        ## Increase brightness manually before displaying 
        # To compensate the grayscale to RGB conversion for tkinter
        self.img1_original = self.adjust_brightness(self.img1_original, brightness_factor)
        self.img2_original = self.adjust_brightness(self.img2_original, brightness_factor)

        self.img1_8bit = img1_8bit
        self.img2_8bit = img2_8bit

        # Canvas size
        self.canvas_width = 600
        self.canvas_height = 400

        # Compute the correct scale factor without stretching
        self.scale_factor = self.get_scale_factor(self.img1_original, self.img2_original)

        # Resize both images using the same scale factor
        self.img1_display = self.resize_image(self.img1_original, self.scale_factor).convert("RGBA")
        self.img2_display = self.resize_image(self.img2_original, self.scale_factor).convert("RGBA")



        # Store the current transformed img1 (for rotation & alpha blending)
        self.img1_transformed = self.img1_display.copy()

        # Convert to Tkinter images
        self.img1_tk = ImageTk.PhotoImage(self.img1_transformed)
        self.img2_tk = ImageTk.PhotoImage(self.img2_display)

        # Create canvas
        self.canvas = tk.Canvas(master, width=self.canvas_width, height=self.canvas_height)
        self.canvas.pack()

        # Display img2_8bit as background (static)
        self.canvas.create_image(0, 0, anchor=tk.NW, image=self.img2_tk)

        # Display img1_8bit (movable)
        self.img1_canvas = self.canvas.create_image(0, 0, anchor=tk.NW, image=self.img1_tk)

        # Alpha slider
        self.alpha_slider = tk.Scale(self.master, from_=0, to=1, resolution=0.01,
                                     orient=tk.HORIZONTAL, label="Alpha", command=self.update_alpha)
        self.alpha_slider.set(0.5)
        self.alpha_slider.pack(fill=tk.X, padx=10, pady=5)

        # Rotation slider (Now Larger & More Precise)
        self.rotation_slider = tk.Scale(self.master, from_=-180, to=180, resolution=0.1, 
                                        orient=tk.HORIZONTAL, label="Rotation", length=500,
                                        command=self.update_rotation)
        self.rotation_slider.set(0)
        self.rotation_slider.pack(fill=tk.X, padx=10, pady=5)

        # Dragging variables
        self.dragging = False
        self.last_x = 0
        self.last_y = 0
        self.img1_x = 0
        self.img1_y = 0

        # Bind mouse events
        self.canvas.bind("<B1-Motion>", self.drag_image)
        self.canvas.bind("<ButtonPress-1>", self.start_drag)
        self.canvas.bind("<ButtonRelease-1>", self.stop_drag)

        # Store alpha & rotation values
        self.alpha = 0.5
        self.angle = 0

        # Store translation values (initially at zero position)
        self.trans_x = 0
        self.trans_y = 0

        # Create Save and Next buttons
        self.save_button = Button(master, text="Save Alignment", command=self.save_alignment)
        self.save_button.pack(fill=tk.X, padx=10, pady=5)

        self.close_button = Button(master, text="Next", command=self.close_window)
        self.close_button.pack(fill=tk.X, padx=10, pady=5)

    def adjust_brightness(self, pil_image, factor):
        """Increase brightness of grayscale image while preserving mode."""
        enhancer = ImageEnhance.Brightness(pil_image)
        return enhancer.enhance(factor)  # Factor > 1 brightens, < 1 darkens

    def get_scale_factor(self, img1, img2):
        """Compute scale factor so that both images fit within the canvas, preserving aspect ratio."""
        w1, h1 = img1.size
        w2, h2 = img2.size

        # Determine the max dimensions of either image
        max_img_width = max(w1, w2)
        max_img_height = max(h1, h2)

        # Calculate scale factors for width and height to fit canvas
        scale_w = self.canvas_width / max_img_width
        scale_h = self.canvas_height / max_img_height

        # Return the smaller of the two to ensure both images fit
        return min(scale_w, scale_h)


    def resize_image(self, img, scale_factor):
        """Resize the image using a uniform scale factor (no distortion)."""
        new_width = int(img.width * scale_factor)
        new_height = int(img.height * scale_factor)
        return img.resize((new_width, new_height), Image.Resampling.LANCZOS)

    def update_alpha(self, event=None):
        """Update transparency of img1."""
        self.alpha = self.alpha_slider.get()  # Store alpha value

        img_with_alpha = self.img1_transformed.copy()  # Apply transparency to transformed img1
        alpha_layer = img_with_alpha.split()[3].point(lambda p: int(p * self.alpha))
        img_with_alpha.putalpha(alpha_layer)

        self.img1_tk = ImageTk.PhotoImage(img_with_alpha)
        self.canvas.itemconfig(self.img1_canvas, image=self.img1_tk)
        
    def update_rotation(self, event=None):
        """Rotate img1 with expand=True and keep it centered properly."""
        self.angle = self.rotation_slider.get()

        # Get current center position of the original image (before rotation)
        orig_center_x = self.img1_display.width // 2
        orig_center_y = self.img1_display.height // 2

        # Rotate with expand=True
        rotated = self.img1_display.rotate(self.angle, resample=Image.Resampling.BICUBIC, expand=True)

        # Calculate offset to preserve center alignment
        new_center_x = rotated.width // 2
        new_center_y = rotated.height // 2

        offset_x = orig_center_x - new_center_x
        offset_y = orig_center_y - new_center_y

        # Apply alpha blending
        alpha_layer = rotated.split()[3].point(lambda p: int(p * self.alpha))
        rotated.putalpha(alpha_layer)

        self.img1_transformed = rotated
        self.img1_tk = ImageTk.PhotoImage(rotated)

        # Update the canvas image
        self.canvas.itemconfig(self.img1_canvas, image=self.img1_tk)

        # Apply offset to keep the center fixed + user's dragging offset
        self.canvas.coords(
            self.img1_canvas,
            self.img1_x + offset_x,
            self.img1_y + offset_y
        )

        # Store offsets in case you want to reference them later
        self.rotation_offset_x = offset_x
        self.rotation_offset_y = offset_y


    def start_drag(self, event):
        self.dragging = True
        self.last_x = event.x
        self.last_y = event.y

    def drag_image(self, event):
        if self.dragging:
            dx = event.x - self.last_x
            dy = event.y - self.last_y

            self.img1_x += dx
            self.img1_y += dy
            self.canvas.move(self.img1_canvas, dx, dy)

            self.last_x = event.x
            self.last_y = event.y

    def stop_drag(self, event):
        self.dragging = False

    def save_alignment(self):
        """Save the transformed image (without blending) to self.manually_aligned_image1."""
        # Get the current transformation values
        self.angle = self.rotation_slider.get()
        self.trans_x = self.img1_x
        self.trans_y = self.img1_y

        # Adjust translation based on the scale factor
        # Scaling the translation values according to the resized images
        self.trans_x /= self.scale_factor
        self.trans_y /= self.scale_factor

        # Original shape
        self.orig_shape = self.img1_8bit.shape
        orig_center = np.array([self.orig_shape[0] / 2, self.orig_shape[1] / 2])
        
        # Apply rotation and translation to the original img1_8bit array (grayscale image)
        # Rotate the original img1_8bit (grayscale)
        rotated_img1 = rotate(self.img1_8bit, self.angle, reshape=True)
        new_shape = rotated_img1.shape
        new_center = np.array([new_shape[0] / 2, new_shape[1] / 2])
        # Shift coordinates to new center
        self.offset = new_center - orig_center

        # Apply offset from rotation with TRUE reshape
        dy, dx = self.offset
        # Apply translation (displacement) to the rotated image (this will shift the image)
        self.manually_aligned_image1 = shift(rotated_img1, (self.trans_y - dy, self.trans_x - dx), mode='nearest')

        # Save the transformation parameters (rotation and translation)
        self.transformation_params = {
            'rotation': self.angle,
            'translation': (self.trans_x, self.trans_y)
        }

        print("Alignment saved.")
        print("Transformation parameters:", self.transformation_params)

    def close_window(self):
        """Close the window."""
        # Explicitly destroy the PhotoImage references to prevent resource leaks
        self.img1_tk = None
        self.img2_tk = None

        # Explicitly delete canvas items (clear the canvas)
        self.canvas.delete(self.img1_canvas)

        # Properly quit and destroy the Tkinter window
        self.master.quit()  # Stop the Tkinter main loop gracefully
        self.master.destroy()  # Destroy the window to free up resources

        # Force garbage collection (after quitting Tkinter and clearing up)
        gc.collect()