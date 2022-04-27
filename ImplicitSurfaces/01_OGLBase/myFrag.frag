#version 330 core

// pipeline-ból bejövõ per-fragment attribútumok
in vec3 vs_out_pos;
in vec3 vs_out_norm;
in vec2 vs_out_tex;
in float vs_out_curv;

out vec4 fs_out_col;

// irány fényforrás: fény iránya
uniform vec3 light_dir = vec3(-1,-1,-1);

// fénytulajdonságok: ambiens, diffúz, ...
uniform vec3 La = vec3(0.1,0.1,0.1);
uniform vec3 Ld = vec3(0.3,0.3,0.3);

uniform sampler2D texImage;
uniform vec3 viewPos;
uniform int showCurv;

uniform vec3 point_light_location = vec3(0,2,1);
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
	

	if(showCurv == 1){
	
		float c = vs_out_curv;
		float a = (1 - c) / 0.25;
		float x = floor(a);
		float y = floor(255 * (a - x));

		float r = 0;
		float g = 0;
		float b = 0;

		if(x == 0){r = 255; g = y; b = 0;}
		if(x == 1){r = 255 - y; g = 255; b = 0;}
		if(x == 2){r = 0; g = 255; b = y;}
		if(x == 3){r = 0; g = 255 - y; b = 255;}
		if(x == 4){r = 0; g = 0; b = 255;}

		if(c == 0){
		
			r = 255;
		
		}

		fs_out_col = vec4(r / 255, g / 255, b / 255, 1);
		
	}
	else {

		fs_out_col = vec4(ambient + diffuse + specular , 1) * texture(texImage, vs_out_tex);

	}



}