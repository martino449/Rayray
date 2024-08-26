import numpy as np
import math
import pygame
from pygame.locals import QUIT
from PIL import Image

# Definizione dei colori
def create_color(r, g, b):
    return np.array([r, g, b])

# Classe per i raggi
class Ray:
    def __init__(self, origin, direction):
        self.origin = np.array(origin)
        self.direction = np.array(direction) / np.linalg.norm(direction)

# Classe per le sfere
class Sphere:
    def __init__(self, center, radius, color):
        self.center = np.array(center)
        self.radius = radius
        self.color = np.array(color)

    def intersect(self, ray):
        oc = ray.origin - self.center
        a = np.dot(ray.direction, ray.direction)
        b = 2.0 * np.dot(oc, ray.direction)
        c = np.dot(oc, oc) - self.radius * self.radius
        discriminant = b * b - 4 * a * c

        if discriminant < 0:
            return False, None
        else:
            t = (-b - math.sqrt(discriminant)) / (2.0 * a)
            if t < 0:
                t = (-b + math.sqrt(discriminant)) / (2.0 * a)
            if t < 0:
                return False, None
            hit_point = ray.origin + t * ray.direction
            normal = (hit_point - self.center) / self.radius
            return True, (hit_point, normal, self.color)

# Classe per i piani
class Plane:
    def __init__(self, point, normal, color):
        self.point = np.array(point)
        self.normal = np.array(normal)
        self.color = np.array(color)

    def intersect(self, ray):
        denom = np.dot(self.normal, ray.direction)
        if abs(denom) > 1e-6:
            t = np.dot(self.point - ray.origin, self.normal) / denom
            if t >= 0:
                hit_point = ray.origin + t * ray.direction
                return True, (hit_point, self.normal, self.color)
        return False, None

# Classe per i cilindri
class Cylinder:
    def __init__(self, center, radius, height, color):
        self.center = np.array(center)
        self.radius = radius
        self.height = height
        self.color = np.array(color)

    def intersect(self, ray):
        co = ray.origin - self.center
        a = ray.direction[0]**2 + ray.direction[2]**2
        b = 2 * (co[0] * ray.direction[0] + co[2] * ray.direction[2])
        c = co[0]**2 + co[2]**2 - self.radius**2
        discriminant = b**2 - 4 * a * c

        if discriminant < 0:
            return False, None

        sqrt_disc = math.sqrt(discriminant)
        t1 = (-b - sqrt_disc) / (2 * a)
        t2 = (-b + sqrt_disc) / (2 * a)

        if t1 > t2:
            t1, t2 = t2, t1

        y1 = ray.origin[1] + t1 * ray.direction[1]
        y2 = ray.origin[1] + t2 * ray.direction[1]

        if y1 < self.center[1] - self.height / 2 or y1 > self.center[1] + self.height / 2:
            if y2 < self.center[1] - self.height / 2 or y2 > self.center[1] + self.height / 2:
                return False, None

        if y1 >= self.center[1] - self.height / 2 and y1 <= self.center[1] + self.height / 2:
            t = t1
        else:
            t = t2

        hit_point = ray.origin + t * ray.direction
        normal = np.array([hit_point[0] - self.center[0], 0, hit_point[2] - self.center[2]])
        normal = normal / np.linalg.norm(normal)

        return True, (hit_point, normal, self.color)

# Classe per la camera
class Camera:
    def __init__(self, position, look_at, up_vector=[0, 1, 0], fov=math.radians(60), aspect_ratio=16/9):
        self.position = np.array(position)
        self.look_at = np.array(look_at)
        self.up_vector = np.array(up_vector)
        self.fov = fov
        self.aspect_ratio = aspect_ratio
        self._update_camera()

    def _update_camera(self):
        self.forward = (self.look_at - self.position) / np.linalg.norm(self.look_at - self.position)
        self.right = np.cross(self.up_vector, self.forward)
        self.right = self.right / np.linalg.norm(self.right)
        self.up = np.cross(self.forward, self.right)
        self.up = self.up / np.linalg.norm(self.up)
        self.half_width = np.tan(self.fov / 2)
        self.half_height = self.half_width / self.aspect_ratio

    def get_ray(self, x, y):
        u = (2 * ((x + 0.5) / width) - 1) * self.half_width
        v = (1 - 2 * ((y + 0.5) / height)) * self.half_height
        direction = self.forward + u * self.right + v * self.up
        return Ray(self.position, direction)

# Funzione per tracciare i raggi
def trace_ray(ray, scene):
    nearest_t = float('inf')
    hit_obj = None
    hit_point = None
    hit_normal = None
    hit_color = None

    for obj in scene:
        hit, data = obj.intersect(ray)
        if hit:
            point, normal, color = data
            t = np.linalg.norm(point - ray.origin)
            if t < nearest_t:
                nearest_t = t
                hit_obj = obj
                hit_point = point
                hit_normal = normal
                hit_color = color

    if hit_obj is None:
        return create_color(0, 0, 0)  # Sfondo nero
    else:
        return hit_color

# Funzione principale per il rendering
def render(camera, scene, width, height):
    image = np.zeros((height, width, 3))
    
    for i in range(height):
        for j in range(width):
            ray = camera.get_ray(j, i)
            pixel_color = trace_ray(ray, scene)
            image[i, j] = np.clip(pixel_color, 0, 1)
        
        # Aggiornamento dell'immagine in tempo reale con pygame
        surface = pygame.surfarray.make_surface((image * 255).astype(np.uint8).swapaxes(0, 1))
        screen.blit(surface, (0, 0))
        pygame.display.flip()
    
    return image

# Impostazioni della scena
rosso = create_color(1, 0, 0)
verde = create_color(0, 1, 0)
blu = create_color(0, 0, 1)
grigio = create_color(0.5, 0.5, 0.5)
giallo = create_color(1, 1, 0)
magenta = create_color(1, 0, 1)

scene = [
    Sphere(center=[2, 0, -6], radius=1, color=verde),  # Sfera verde
    Sphere(center=[-2, 0, -6], radius=1, color=blu),  # Sfera blu
    Plane(point=[0, -1, 0], normal=[0, 1, 0], color=grigio),  # Piano grigio
    Cylinder(center=[-4, 0, -6], radius=1, height=2, color=giallo),  # Cilindro giallo
]

# Impostazioni della camera e della finestra di pygame
pygame.init()
width, height = 600, 400
screen = pygame.display.set_mode((width, height))

# Rendering con diverse angolazioni della camera
scena = 1
while scena < 20:
    camera1 = Camera(
        position=[1, scena, 5],  # Position the camera behind the objects
        look_at=[0 , 1, -5],  # Look towards the center of the scene
        aspect_ratio = width / height
    )
    scena += 0.5
    
    # Rendering e visualizzazione in tempo reale
    render(camera1, scene, width, height)
    
    # Salvataggio dell'immagine
    pygame.image.save(screen, f'render_camera_t3-{scena}.png')
    
    # Gestione degli eventi
    for event in pygame.event.get():
        if event.type == QUIT:
            pygame.quit()
            exit()

pygame.quit()
