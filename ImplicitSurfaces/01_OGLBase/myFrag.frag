#version 330 core

// pipeline-ból bejövõ per-fragment attribútumok
in vec3 vs_out_pos;
in vec3 vs_out_norm;
in vec2 vs_out_tex;

out vec4 fs_out_col;

// irány fényforrás: fény iránya
uniform vec3 light_dir = vec3(-1,-1,-1);

// fénytulajdonságok: ambiens, diffúz, ...
uniform vec3 La = vec3(0.1,0.1,0.1);
uniform vec3 Ld = vec3(0.3,0.3,0.3);

uniform sampler2D texImage;

uniform vec3 viewPos;

uniform vec3 point_light_location = vec3(0,2,0);
uniform vec3 point_light_color = vec3(0.4,0.4,0.4);


void main()
{

	vec3 ambient = La;
	vec3 normal = normalize(vs_out_norm);
	
	vec3 point_light_direction = normalize(vs_out_pos - point_light_location);

	vec3 to_light = normalize(-light_dir);
	vec3 to_point_light = point_light_direction * -1;
	

	float cosa = clamp(dot(normal, to_light), 0, 1);
	float cosa_point = clamp(dot(normal, to_point_light), 0, 1);


	vec3 diffuse = cosa*Ld;
	vec3 diffuse_point = cosa_point * point_light_color;

	float specularStrength = 1;

	vec3 viewDir = normalize(viewPos - vs_out_pos);

	vec3 reflectDir = reflect(light_dir, normal);  
	vec3 reflectDir_point = reflect(-to_point_light, normal);

	vec3 h = normalize(viewDir + reflectDir);
	vec3 h_point = normalize(viewDir + reflectDir_point);

	float spec = pow(max(dot(viewDir, h), 0.0), 48);
	float spec_point = pow(max(dot(viewDir, h_point), 0.0), 32);
	vec3 lightColor = vec3(0.3,0.3,0.3);

	vec3 specular = specularStrength * spec * lightColor; 
	vec3 specular_point = specularStrength * spec_point * point_light_color * 0;

	fs_out_col = vec4(ambient + diffuse + diffuse_point + specular + specular_point, 1) * texture(texImage, vs_out_tex);
}