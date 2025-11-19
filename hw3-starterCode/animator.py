import subprocess
import math
from moviepy import ImageSequenceClip
import os

# the general idea is to create an infinity animation using
# Lemniscate of Gerono and then have this python script run
# the cpp program for me automatically
def genImages():
   i = 1
   a = 17
   step_size = (2 * math.pi)/100
   t1 = (math.pi/2)
   t2 = t1 + (step_size * 12)
   t3 = t2 + (step_size * 12)
   while i <= 100: # 100 frames
      t1 += step_size
      t2 += step_size
      t3 += step_size

      # final x, y positions
      x1 = a * math.sin(t1)
      y1 = x1 * math.cos(t1)
      x2 = a * math.sin(t2)
      y2 = x2 * math.cos(t2)
      x3 = a * math.sin(t3)
      y3 = x2 * math.cos(t3)

      with open(f'animation/animation_{i}.scene', 'w') as f:
         f.write(f"""4
         amb: 0.100000 0.100000 0.100000
         sphere
         pos: {x1} {y1} -40.000000
         rad: 3.000000
         dif: 0.0 0.0 1.0
         spe: 0.500000 0.500000 0.500000
         shi: 10.000000
         sphere
         pos: {x2} {y2} -40.000000
         rad: 3.000000
         dif: 1.0 0.0 0.0
         spe: 0.500000 0.500000 0.500000
         shi: 10.000000
         sphere
         pos: {x3} {y3} -40.000000
         rad: 3.000000
         dif: 0.0 1.0 0.0
         spe: 0.500000 0.500000 0.500000
         shi: 10.000000
         light
         pos: 0.000000 10.000000 -35.000000
         col: 1.000000 1.000000 1.000000
         """)
         f.close()
      result = subprocess.run(['./hw3', f'animation/animation_{i}.scene', f'animation/frame_{i}.jpg'])
      i = i + 1

def images_to_mp4(image_folder, output_file, fps):
    images = [img for img in os.listdir(image_folder)
      if img.lower().endswith((".png", ".jpg", ".jpeg"))]
    images.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))
    
    images = [os.path.join(image_folder, img) for img in images]
    
    clip = ImageSequenceClip(images, fps=fps)
    clip.write_videofile(output_file, codec="libx264")

if __name__ == "__main__":
   genImages()
   images_to_mp4("animation", "final_animation.mp4", 25)
