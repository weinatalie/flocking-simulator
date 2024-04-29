#version 330

uniform vec4 u_color;
uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;
uniform bool skybox;
uniform bool predator;
uniform samplerCube u_texture_cubemap;
uniform sampler2D u_texture_1;
uniform sampler2D u_texture_2;

in vec4 v_position;
in vec4 v_normal;
in vec2 v_uv;

out vec4 out_color;

void main() {
    float p = 20.0;
     float ka = 0.25;
     float Ia = 2.0;
     float ks = 0.4;
     float kd = 0.6;

     float r = distance(u_light_pos, v_position.xyz);
     vec3 n = normalize(v_normal.xyz);
     vec3 l = normalize(u_light_pos - v_position.xyz);
     vec3 h_v = normalize(l + normalize(u_cam_pos - v_position.xyz));
     float ambient = ka * Ia;
     vec3 diffuse = kd * (u_light_intensity/(r*r)) * max(0, dot(n, l));
     vec3 specular = ks * (u_light_intensity/(r*r)) * pow(max(0, dot(n, h_v)), p);
   out_color = vec4(ambient + diffuse + specular, 0);
    out_color.a = 1;

    if (skybox) {
    mediump vec3 w_o = normalize(v_position.xyz / v_position.w).xyz - u_cam_pos.xyz;
    out_color.xyz = texture(u_texture_cubemap, w_o).xyz;
  } else if (predator) {
    vec3 tex_color = texture(u_texture_2, v_uv).xyz;
    vec3 shading = ambient + diffuse + specular;
    out_color = vec4(shading * tex_color, 1.0);
  } else {
    vec3 tex_color = texture(u_texture_1, v_uv).xyz;
    vec3 shading = ambient + diffuse + specular;
    out_color = vec4(shading * tex_color, 1.0);
  }
}