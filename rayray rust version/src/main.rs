use image::{ImageBuffer, Rgb};
use std::f64::consts::PI;

// Funzione per creare colori
fn create_color(r: f64, g: f64, b: f64) -> [f64; 3] {
    [r, g, b]
}

// Struttura per i raggi
struct Ray {
    origin: [f64; 3],
    direction: [f64; 3],
}

impl Ray {
    fn new(origin: [f64; 3], direction: [f64; 3]) -> Self {
        let norm = (direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2])
            .sqrt();
        Self {
            origin,
            direction: [direction[0] / norm, direction[1] / norm, direction[2] / norm],
        }
    }
}

// Struttura per le sfere
struct Sphere {
    center: [f64; 3],
    radius: f64,
    color: [f64; 3],
}

impl Sphere {
    fn new(center: [f64; 3], radius: f64, color: [f64; 3]) -> Self {
        Self {
            center,
            radius,
            color,
        }
    }

    fn intersect(&self, ray: &Ray) -> Option<([f64; 3], [f64; 3], [f64; 3])> {
        let oc = [
            ray.origin[0] - self.center[0],
            ray.origin[1] - self.center[1],
            ray.origin[2] - self.center[2],
        ];

        let a = ray.direction.iter().map(|&x| x * x).sum::<f64>();
        let b = 2.0 * oc.iter().zip(ray.direction.iter()).map(|(&o, &d)| o * d).sum::<f64>();
        let c = oc.iter().map(|&o| o * o).sum::<f64>() - self.radius * self.radius;
        let discriminant = b * b - 4.0 * a * c;

        if discriminant < 0.0 {
            return None;
        }

        let t1 = (-b - discriminant.sqrt()) / (2.0 * a);
        let t2 = (-b + discriminant.sqrt()) / (2.0 * a);

        let t = if t1 < 0.0 { t2 } else { t1 };

        if t < 0.0 {
            return None;
        }

        let hit_point = [
            ray.origin[0] + t * ray.direction[0],
            ray.origin[1] + t * ray.direction[1],
            ray.origin[2] + t * ray.direction[2],
        ];

        let normal = [
            (hit_point[0] - self.center[0]) / self.radius,
            (hit_point[1] - self.center[1]) / self.radius,
            (hit_point[2] - self.center[2]) / self.radius,
        ];

        Some((hit_point, normal, self.color))
    }
}

// Struttura per i piani
struct Plane {
    point: [f64; 3],
    normal: [f64; 3],
    color: [f64; 3],
}

impl Plane {
    fn new(point: [f64; 3], normal: [f64; 3], color: [f64; 3]) -> Self {
        Self { point, normal, color }
    }

    fn intersect(&self, ray: &Ray) -> Option<([f64; 3], [f64; 3], [f64; 3])> {
        let denom = self.normal.iter().zip(ray.direction.iter()).map(|(&n, &d)| n * d).sum::<f64>();

        if denom.abs() > 1e-6 {
            let t = self.point.iter()
                .zip(ray.origin.iter())
                .zip(self.normal.iter())
                .map(|((&p, &o), &n)| (p - o) * n)
                .sum::<f64>() / denom;

            if t >= 0.0 {
                let hit_point = [
                    ray.origin[0] + t * ray.direction[0],
                    ray.origin[1] + t * ray.direction[1],
                    ray.origin[2] + t * ray.direction[2],
                ];

                return Some((hit_point, self.normal, self.color));
            }
        }

        None
    }
}

// Struttura per i cilindri
struct Cylinder {
    center: [f64; 3],
    radius: f64,
    height: f64,
    color: [f64; 3],
}

impl Cylinder {
    fn new(center: [f64; 3], radius: f64, height: f64, color: [f64; 3]) -> Self {
        Self {
            center,
            radius,
            height,
            color,
        }
    }

    fn intersect(&self, ray: &Ray) -> Option<([f64; 3], [f64; 3], [f64; 3])> {
        let co = [
            ray.origin[0] - self.center[0],
            ray.origin[1] - self.center[1],
            ray.origin[2] - self.center[2],
        ];

        let a = ray.direction[0] * ray.direction[0] + ray.direction[2] * ray.direction[2];
        let b = 2.0 * (co[0] * ray.direction[0] + co[2] * ray.direction[2]);
        let c = co[0] * co[0] + co[2] * co[2] - self.radius * self.radius;

        let discriminant = b * b - 4.0 * a * c;

        if discriminant < 0.0 {
            return None;
        }

        let sqrt_disc = discriminant.sqrt();
        let mut t1 = (-b - sqrt_disc) / (2.0 * a);
        let mut t2 = (-b + sqrt_disc) / (2.0 * a);

        if t1 > t2 {
            std::mem::swap(&mut t1, &mut t2);
        }

        let y1 = ray.origin[1] + t1 * ray.direction[1];
        let y2 = ray.origin[1] + t2 * ray.direction[1];

        if (y1 < self.center[1] - self.height / 2.0) || (y1 > self.center[1] + self.height / 2.0) {
            if (y2 < self.center[1] - self.height / 2.0) || (y2 > self.center[1] + self.height / 2.0) {
                return None;
            }
        }

        let t = if y1 >= self.center[1] - self.height / 2.0 && y1 <= self.center[1] + self.height / 2.0 {
            t1
        } else {
            t2
        };

        let hit_point = [
            ray.origin[0] + t * ray.direction[0],
            ray.origin[1] + t * ray.direction[1],
            ray.origin[2] + t * ray.direction[2],
        ];

        let normal = [
            hit_point[0] - self.center[0],
            0.0,
            hit_point[2] - self.center[2],
        ];
        let normal_norm = (normal[0] * normal[0] + normal[2] * normal[2]).sqrt();

        Some((
            hit_point,
            [normal[0] / normal_norm, 0.0, normal[2] / normal_norm],
            self.color,
        ))
    }
}

// Struttura per la camera
struct Camera {
    position: [f64; 3],
    look_at: [f64; 3],
    up_vector: [f64; 3],
    fov: f64,
    aspect_ratio: f64,
    forward: [f64; 3],
    right: [f64; 3],
    up: [f64; 3],
    half_width: f64,
    half_height: f64,
}

impl Camera {
    fn new(position: [f64; 3], look_at: [f64; 3], up_vector: [f64; 3], fov: f64, aspect_ratio: f64) -> Self {
        let forward = [
            look_at[0] - position[0],
            look_at[1] - position[1],
            look_at[2] - position[2],
        ];
        let norm = (forward[0] * forward[0] + forward[1] * forward[1] + forward[2] * forward[2]).sqrt();
        let forward = [forward[0] / norm, forward[1] / norm, forward[2] / norm];

        let right = [
            up_vector[1] * forward[2] - up_vector[2] * forward[1],
            up_vector[2] * forward[0] - up_vector[0] * forward[2],
            up_vector[0] * forward[1] - up_vector[1] * forward[0],
        ];
        let norm = (right[0] * right[0] + right[1] * right[1] + right[2] * right[2]).sqrt();
        let right = [right[0] / norm, right[1] / norm, right[2] / norm];

        let up = [
            forward[1] * right[2] - forward[2] * right[1],
            forward[2] * right[0] - forward[0] * right[2],
            forward[0] * right[1] - forward[1] * right[0],
        ];
        let norm = (up[0] * up[0] + up[1] * up[1] + up[2] * up[2]).sqrt();
        let up = [up[0] / norm, up[1] / norm, up[2] / norm];

        let half_width = (fov / 2.0).tan();
        let half_height = half_width / aspect_ratio;

        Self {
            position,
            look_at,
            up_vector,
            fov,
            aspect_ratio,
            forward,
            right,
            up,
            half_width,
            half_height,
        }
    }

    fn get_ray(&self, x: f64, y: f64, width: f64, height: f64) -> Ray {
        let u = (2.0 * ((x + 0.5) / width) - 1.0) * self.half_width;
        let v = (1.0 - 2.0 * ((y + 0.5) / height)) * self.half_height;

        let direction = [
            self.forward[0] + u * self.right[0] + v * self.up[0],
            self.forward[1] + u * self.right[1] + v * self.up[1],
            self.forward[2] + u * self.right[2] + v * self.up[2],
        ];

        Ray::new(self.position, direction)
    }
}

// Funzione per tracciare i raggi
fn trace_ray(ray: &Ray, scene: &[Box<dyn Intersectable>]) -> [f64; 3] {
    let mut nearest_t = f64::INFINITY;
    let mut hit_color = [0.0, 0.0, 0.0]; // Nero per il background

    for obj in scene {
        if let Some((_, _, color)) = obj.intersect(ray) {
            let t = (ray.origin.iter().zip(color.iter()).map(|(&o, &c)| (o - c).abs()).sum::<f64>()).sqrt();
            if t < nearest_t {
                nearest_t = t;
                hit_color = color;
            }
        }
    }

    hit_color
}

// Funzione principale per il rendering
fn render(camera: &Camera, scene: &[Box<dyn Intersectable>], width: u32, height: u32) -> ImageBuffer<Rgb<u8>, Vec<u8>> {
    let mut img = ImageBuffer::new(width, height);

    for y in 0..height {
        for x in 0..width {
            let ray = camera.get_ray(x as f64, y as f64, width as f64, height as f64);
            let color = trace_ray(&ray, scene);
            let pixel_color = [
                (color[0].min(1.0) * 255.0) as u8,
                (color[1].min(1.0) * 255.0) as u8,
                (color[2].min(1.0) * 255.0) as u8,
            ];
            img.put_pixel(x, y, Rgb(pixel_color));
        }
    }

    img
}

// Trait per oggetti che possono essere intersecati
trait Intersectable {
    fn intersect(&self, ray: &Ray) -> Option<([f64; 3], [f64; 3], [f64; 3])>;
}

// Implementazione del trait per Sphere
impl Intersectable for Sphere {
    fn intersect(&self, ray: &Ray) -> Option<([f64; 3], [f64; 3], [f64; 3])> {
        Sphere::intersect(self, ray)
    }
}

// Implementazione del trait per Plane
impl Intersectable for Plane {
    fn intersect(&self, ray: &Ray) -> Option<([f64; 3], [f64; 3], [f64; 3])> {
        Plane::intersect(self, ray)
    }
}

// Implementazione del trait per Cylinder
impl Intersectable for Cylinder {
    fn intersect(&self, ray: &Ray) -> Option<([f64; 3], [f64; 3], [f64; 3])> {
        Cylinder::intersect(self, ray)
    }
}

fn main() {
    let verde = create_color(0.0, 1.0, 0.0);
    let blu = create_color(0.0, 0.0, 1.0);
    let grigio = create_color(0.5, 0.5, 0.5);
    let giallo = create_color(1.0, 1.0, 0.0);

    let scene: Vec<Box<dyn Intersectable>> = vec![
        Box::new(Sphere::new([2.0, 0.0, -6.0], 1.0, verde)),
        Box::new(Sphere::new([-2.0, 0.0, -6.0], 1.0, blu)),
        Box::new(Plane::new([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], grigio)),
        Box::new(Cylinder::new([-4.0, 0.0, -6.0], 1.0, 2.0, giallo)),
    ];

    let width = 600;
    let height = 400;

    let mut scena = 1.0;

    while scena < 50.0 {
        let camera = Camera::new(
            [1.0, scena, 5.0],
            [0.0, 1.0, -5.0],
            [0.0, 1.0, 0.0],
            60.0 * PI / 180.0,
            width as f64 / height as f64,
        );

        let img = render(&camera, &scene, width, height);
        img.save(format!("render_camera_t3-{:.1}.png", scena)).unwrap();

        scena += 0.5;
    }
}
